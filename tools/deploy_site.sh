#!/bin/sh
# Build the static site (web UI + index.json + SVG) and publish it to the gh-pages branch.
set -e
cd "$(dirname "$0")/.."
.venv/bin/python -m sdb site
WT=$(mktemp -d)
git fetch -q origin gh-pages 2>/dev/null || true
if git rev-parse -q --verify origin/gh-pages >/dev/null; then
  git worktree add -q -B gh-pages "$WT" origin/gh-pages
else
  git worktree add -q --detach "$WT"
  (cd "$WT" && git checkout -q --orphan gh-pages && git rm -rqf .)
fi
rsync -a --delete --exclude .git site/ "$WT/"
cd "$WT"
git add -A
git commit -q -m "site: $(date '+%Y-%m-%d %H:%M')" || echo "沒有變更"
git push -q origin gh-pages
cd - >/dev/null
git worktree remove --force "$WT"
echo "已推到 gh-pages"
