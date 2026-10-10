"""Development server.

Sends Cache-Control: no-cache, and stamps the scripts and style.css with their
modification times as they are referenced from index.html.  The header alone is
not enough: an embedded webview such as VS Code's Simple Browser keeps its own
cache and will happily run yesterday's script against today's markup, which
looks exactly like a feature that has stopped working.
"""
import http.server
import os
import re
import sys


class H(http.server.SimpleHTTPRequestHandler):

    def end_headers(self):
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        http.server.SimpleHTTPRequestHandler.end_headers(self)

    def do_GET(self):
        path = self.path.split("?")[0]
        if path in ("/", "/index.html"):
            return self._index()
        return http.server.SimpleHTTPRequestHandler.do_GET(self)

    def _index(self):
        try:
            html = open("index.html", encoding="utf-8").read()
        except OSError:
            self.send_error(404)
            return

        def stamp(m):
            name = m.group(1)
            try:
                t = int(os.path.getmtime(name))
            except OSError:
                t = 0
            return '%s?t=%d' % (name, t)

        html = re.sub(r'\b(app\.js|physics\.js|style\.css)\b', stamp, html)
        body = html.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    print("http://localhost:%d/" % port)
    http.server.HTTPServer(("", port), H).serve_forever()
