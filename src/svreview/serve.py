"""serve.py - a local web server for a site directory: range requests (the BAM reader asks for the bytes it
needs) and verdicts written to disk as they are given.

  GET  <anything>                 the file, whole or the requested byte range (206, Accept-Ranges)
  GET  <dir>/verdicts.tsv         the verdicts so far (404 before the first)
  POST <dir>/verdicts.tsv         one or more verdict rows (the TSV body, with or without its header) appended
                                  to that file; the response is the file's row count

A static host (GitHub Pages) has no POST: the viewer then keeps the verdicts in the browser and offers the
TSV as a download. Python's own http.server has no range support, which is why this exists.
"""

import email.utils
import http.server
import os
import posixpath
import re
import sys
import urllib.parse

from . import contract

RANGE_RE = re.compile(r"^bytes=(\d*)-(\d*)$")


class Handler(http.server.SimpleHTTPRequestHandler):
    server_version = "svreview/0.1"
    extensions_map = dict(http.server.SimpleHTTPRequestHandler.extensions_map)
    extensions_map.update({".bam": "application/octet-stream", ".bai": "application/octet-stream",
                           ".json": "application/json", ".tsv": "text/tab-separated-values", ".mjs": "text/javascript",
                           ".js": "text/javascript"})

    def log_message(self, fmt, *args):
        if self.server.quiet:
            return
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-cache")
        http.server.SimpleHTTPRequestHandler.end_headers(self)

    def do_GET(self):
        path = self.translate_path(self.path)
        rng = self.headers.get("Range")
        if rng and os.path.isfile(path):
            return self._send_range(path, rng)
        return http.server.SimpleHTTPRequestHandler.do_GET(self)

    def _send_range(self, path, rng):
        m = RANGE_RE.match(rng.strip())
        size = os.path.getsize(path)
        if not m:
            self.send_error(416, "bad range")
            return
        first, last = m.group(1), m.group(2)
        if first == "" and last == "":
            self.send_error(416, "bad range")
            return
        if first == "":
            n = int(last)
            start, end = max(0, size - n), size - 1
        else:
            start = int(first)
            end = min(int(last), size - 1) if last else size - 1
        if start >= size or start > end:
            self.send_response(416)
            self.send_header("Content-Range", "bytes */%d" % size)
            self.end_headers()
            return
        self.send_response(206)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Range", "bytes %d-%d/%d" % (start, end, size))
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Last-Modified", email.utils.formatdate(os.path.getmtime(path), usegmt=True))
        self.end_headers()
        with open(path, "rb") as fh:
            fh.seek(start)
            remaining = end - start + 1
            while remaining > 0:
                chunk = fh.read(min(65536, remaining))
                if not chunk:
                    break
                self.wfile.write(chunk)
                remaining -= len(chunk)

    def do_POST(self):
        path = self.translate_path(self.path)
        if os.path.basename(path) != "verdicts.tsv":
            self.send_error(404, "only verdicts.tsv takes a POST")
            return
        root = os.path.abspath(self.directory)
        if not os.path.abspath(path).startswith(root + os.sep):
            self.send_error(403, "outside the served directory")
            return
        length = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(length).decode("utf-8") if length else ""
        rows = []
        for line in body.splitlines():
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "candidate" and "verdict" in f:
                continue                      # a header line
            if len(f) < len(contract.VERDICT_COLUMNS):
                f += [""] * (len(contract.VERDICT_COLUMNS) - len(f))
            row = dict(zip(contract.VERDICT_COLUMNS, f))
            if row["verdict"] and row["verdict"] not in contract.VERDICTS:
                self.send_error(400, "verdict %r is not one of %s" % (row["verdict"], ", ".join(contract.VERDICTS)))
                return
            if contract.to_int(row["start"]) is None or contract.to_int(row["end"]) is None:
                self.send_error(400, "start/end are not integers")
                return
            rows.append(row)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        contract.write_verdicts(path, rows, append=True)
        n = len(contract.read_verdicts(path))
        payload = ('{"appended": %d, "verdicts": %d, "path": "%s"}' % (len(rows), n, os.path.basename(path))).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def translate_path(self, path):
        # SimpleHTTPRequestHandler's, rooted at the served directory (no query, no fragment)
        path = path.split("?", 1)[0].split("#", 1)[0]
        path = posixpath.normpath(urllib.parse.unquote(path))
        words = [w for w in path.split("/") if w and w not in (os.curdir, os.pardir)]
        out = self.directory
        for w in words:
            out = os.path.join(out, w)
        return out


class Server(http.server.ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def __init__(self, directory, port=8765, host="127.0.0.1", quiet=False):
        self.directory = os.path.abspath(directory)
        self.quiet = quiet
        handler = lambda *a, **kw: Handler(*a, directory=self.directory, **kw)   # noqa: E731
        http.server.ThreadingHTTPServer.__init__(self, (host, port), handler)


def serve(directory, port=8765, host="127.0.0.1", quiet=False, open_browser=False, path="index.html"):
    """Serve *directory* until interrupted."""
    srv = Server(directory, port, host, quiet)
    url = "http://%s:%d/%s" % (host, srv.server_address[1], path)
    print("svreview: serving %s at %s (Ctrl-C stops)" % (srv.directory, url))
    if open_browser:
        import webbrowser
        webbrowser.open(url)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
