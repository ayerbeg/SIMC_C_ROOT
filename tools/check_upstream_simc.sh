#!/usr/bin/env bash
# tools/check_upstream_simc.sh — is the official Fortran SIMC ahead of what SIMC++ follows?
#
# SIMC++ records in UPSTREAM_SIMC the commit of JeffersonLab/simc_gfortran whose behaviour it reproduces
# (--fortran-compat target). This script lists the upstream commits made after it, with the files they touch,
# so a maintainer can decide what to port. Needs either the GitHub CLI (gh) or git with network access.
#
# usage: tools/check_upstream_simc.sh [branch]      (default branch: master)
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
base=$(grep -E '^commit' "$here/UPSTREAM_SIMC" | awk '{print $2}')
repo=JeffersonLab/simc_gfortran
branch=${1:-master}
if command -v gh >/dev/null 2>&1; then
  out=$(gh api "repos/$repo/compare/$base...$branch" --jq '"ahead_by \(.ahead_by)", (.commits[] | "\(.sha[0:7]) \(.commit.author.date[0:10]) \(.commit.message | split("\n")[0])"), "--- files", (.files[]? | "\(.filename) +\(.additions) -\(.deletions)")')
else
  tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
  git clone -q --filter=blob:none --no-checkout "https://github.com/$repo.git" "$tmp/simc"
  out=$(cd "$tmp/simc" && echo "ahead_by $(git rev-list --count "$base..origin/$branch")" && git log --format='%h %ad %s' --date=short "$base..origin/$branch" && echo "--- files" && git diff --stat "$base" "origin/$branch")
fi
echo "SIMC++ follows $repo@$base"
echo "$out"
n=$(echo "$out" | awk '/^ahead_by/{print $2}')
if [ "${n:-0}" = 0 ]; then echo "Up to date."; else echo "Upstream has $n newer commit(s): review them (docs/translation, 'Following upstream SIMC')."; fi
