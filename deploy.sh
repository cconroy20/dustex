#!/bin/bash
# Publish web/ to the gh-pages branch as a SINGLE commit, as for spex,
# ismex and ionbal: the branch is rebuilt from scratch and force-pushed on every deploy.
#
#   ./deploy.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BRANCH=gh-pages
cd "$HERE"
git rev-parse --git-dir > /dev/null 2>&1 || { echo "not a git repository"; exit 1; }
git remote get-url origin > /dev/null 2>&1 || { echo "no 'origin' remote"; exit 1; }
[ -f web/index.html ] || { echo "no web/index.html"; exit 1; }
[ -f web/.nojekyll ] || { echo "web/.nojekyll missing"; exit 1; }

# Build on a throwaway branch and push it to gh-pages by refspec.
TMP=$(mktemp -d)
SCRATCH="deploy-$$"
git worktree add --detach -q "$TMP"
(
  cd "$TMP"
  git checkout -q --orphan "$SCRATCH"
  git rm -rqf . 2>/dev/null || true
  cp -R "$HERE"/web/. .
  rm -f serve.py
  # Pages caches index.html and each script independently for ten minutes, so
  # new HTML could run against an old script.  Stamping the commit onto every
  # asset URL makes new HTML always ask for URLs it has never seen.
  V=$(git -C "$HERE" rev-parse --short HEAD)
  for f in app.js physics.js; do
    sed -i '' -e "s|src=\"$f\"|src=\"$f?v=$V\"|" index.html
    grep -q "$f?v=" index.html || { echo "asset stamping failed for $f"; exit 1; }
  done
  sed -i '' -e "s|href=\"style.css\"|href=\"style.css?v=$V\"|" index.html
  grep -q "style.css?v=" index.html || { echo "asset stamping failed for style.css"; exit 1; }
  git add -A
  git commit -q -m "site $(date -u '+%Y-%m-%d %H:%M UTC')"
  git push -f -q origin "HEAD:$BRANCH"
)
git worktree remove --force "$TMP"
git branch -qD "$SCRATCH" 2>/dev/null || true
echo "pushed $BRANCH  ($(du -sh "$HERE"/web | cut -f1))"
