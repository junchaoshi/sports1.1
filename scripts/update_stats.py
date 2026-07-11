#!/usr/bin/env python3
"""Update README traffic statistics from GitHub's 14-day Traffic API."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import re
import sys
import urllib.request
from urllib.parse import quote


REPO = "junchaoshi/sports1.1"
STATS_START = "<!-- stats:start -->"
STATS_END = "<!-- stats:end -->"


def utc_now() -> dt.datetime:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0)


def api_json(path: str, token: str) -> dict:
    req = urllib.request.Request(
        f"https://api.github.com{path}",
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "sports1.1-stats-updater",
        },
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.load(response)


def fetch_traffic(repo: str, metric: str, token: str) -> list[dict]:
    payload = api_json(f"/repos/{repo}/traffic/{metric}?per=day", token)
    return payload[metric]


def load_stats(path: pathlib.Path) -> dict:
    if not path.exists():
        return {"repository": REPO, "daily": {}}
    return json.loads(path.read_text(encoding="utf-8"))


def merge_metric(stats: dict, metric: str, rows: list[dict]) -> None:
    daily = stats.setdefault("daily", {})
    for row in rows:
        day = row["timestamp"][:10]
        entry = daily.setdefault(day, {})
        entry[metric] = {"count": int(row["count"]), "uniques": int(row["uniques"])}


def count_for(day_entry: dict, metric: str) -> int:
    return int(day_entry.get(metric, {}).get("count", 0))


def summarize(stats: dict, today: dt.date) -> dict:
    month_prefix = today.strftime("%Y-%m")
    year_prefix = today.strftime("%Y")
    summary = {
        "month": {"views": 0, "clones": 0},
        "year": {"views": 0, "clones": 0},
        "total": {"views": 0, "clones": 0},
    }
    for day, entry in stats.get("daily", {}).items():
        for metric in ("views", "clones"):
            count = count_for(entry, metric)
            summary["total"][metric] += count
            if day.startswith(year_prefix):
                summary["year"][metric] += count
            if day.startswith(month_prefix):
                summary["month"][metric] += count
    return summary


def comma(value: int) -> str:
    return f"{value:,}"


def badge(label: str, value: int, color: str) -> str:
    label_part = quote(label, safe="")
    value_part = quote(comma(value), safe="")
    color_part = quote(color, safe="")
    return f"![{label}: {comma(value)}](https://img.shields.io/badge/{label_part}-{value_part}-{color_part}?style=flat-square)"


def render_block(stats: dict) -> str:
    summary = stats["summary"]
    updated = stats["updated_at_utc"][:10]
    return "\n".join(
        [
            STATS_START,
            "| Metric | This month | This year | Total |",
            "|---|---|---|---|",
            f"| Repository views | {badge('month', summary['month']['views'], 'blue')} | {badge('year', summary['year']['views'], 'blue')} | {badge('total', summary['total']['views'], 'blue')} |",
            f"| Git clones (download proxy) | {badge('month', summary['month']['clones'], 'brightgreen')} | {badge('year', summary['year']['clones'], 'brightgreen')} | {badge('total', summary['total']['clones'], 'brightgreen')} |",
            "",
            f"Last updated: {updated} UTC. Git clones are used as the download proxy because GitHub does not report `master.zip` download counts.",
            STATS_END,
        ]
    )


def update_readme(text: str, block: str) -> str:
    marked = re.compile(
        rf"{re.escape(STATS_START)}.*?{re.escape(STATS_END)}",
        re.DOTALL,
    )
    if marked.search(text):
        return marked.sub(block, text)

    heading = re.search(r"(?m)^## Software statistics <a id='statistics'></a>\n", text)
    if not heading:
        return text.rstrip() + "\n\n## Software statistics <a id='statistics'></a>\n" + block + "\n"

    next_heading = re.search(r"(?m)^## ", text[heading.end() :])
    end = heading.end() + next_heading.start() if next_heading else len(text)
    return text[: heading.end()] + block + "\n\n" + text[end:].lstrip("\n")


def write_json(path: pathlib.Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def update(repo: str, stats_path: pathlib.Path, readme_path: pathlib.Path, token: str) -> None:
    now = utc_now()
    stats = load_stats(stats_path)
    stats["repository"] = repo
    merge_metric(stats, "views", fetch_traffic(repo, "views", token))
    merge_metric(stats, "clones", fetch_traffic(repo, "clones", token))
    stats["updated_at_utc"] = now.isoformat().replace("+00:00", "Z")
    stats["summary"] = summarize(stats, now.date())
    write_json(stats_path, stats)
    readme_path.write_text(update_readme(readme_path.read_text(encoding="utf-8"), render_block(stats)), encoding="utf-8")


def self_test() -> None:
    stats = {"daily": {"2025-12-31": {"views": {"count": 5, "uniques": 4}}}}
    merge_metric(
        stats,
        "views",
        [
            {"timestamp": "2026-07-01T00:00:00Z", "count": 10, "uniques": 8},
            {"timestamp": "2026-07-02T00:00:00Z", "count": 20, "uniques": 9},
        ],
    )
    merge_metric(stats, "clones", [{"timestamp": "2026-07-02T00:00:00Z", "count": 3, "uniques": 2}])
    summary = summarize(stats, dt.date(2026, 7, 9))
    assert summary == {
        "month": {"views": 30, "clones": 3},
        "year": {"views": 30, "clones": 3},
        "total": {"views": 35, "clones": 3},
    }
    stats["summary"] = summary
    stats["updated_at_utc"] = "2026-07-09T00:00:00Z"
    block = render_block(stats)
    readme = "Intro\n\n## Software statistics <a id='statistics'></a>\nold badge\n\n## Disclaimer <a id='disclaimer'></a>\ntext\n"
    updated = update_readme(readme, block)
    assert "old badge" not in updated
    assert "img.shields.io/badge/month-3-brightgreen" in updated
    assert "## Disclaimer <a id='disclaimer'></a>" in updated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", REPO))
    parser.add_argument("--stats-file", default="metrics/stats.json")
    parser.add_argument("--readme", default="README.md")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return 0

    token = os.environ.get("TRAFFIC_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        print("TRAFFIC_TOKEN or GITHUB_TOKEN is required for GitHub Traffic API access.", file=sys.stderr)
        return 2

    update(args.repo, pathlib.Path(args.stats_file), pathlib.Path(args.readme), token)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
