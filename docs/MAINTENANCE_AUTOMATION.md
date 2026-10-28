# OpenAI Codex maintenance every two hours

The repository includes
[`.github/workflows/codex-maintenance.yml`](../.github/workflows/codex-maintenance.yml),
using the official [Codex GitHub Action](https://learn.chatgpt.com/docs/github-action).
It proposes maintenance through draft pull requests; it never merges changes.

## Activate

1. Review and merge the workflow and companion files into the default branch of
   `pha4ge/CholeraFlow`. The workflow deliberately does not run from other branches
   or forks; update that repository guard explicitly when maintaining a fork.
2. Add repository Actions secret `OPENAI_API_KEY` containing a project API key
   permitted to use Codex. Enter it in GitHub settings or with `gh secret set`;
   never place it in source files, prompts, test output or chat.
3. Enable GitHub Actions and allow Actions to create pull requests under
   Settings → Actions → General. Keep protected-branch/human review requirements.
4. Optionally set repository variable `CODEX_MODEL` to an available model ID.
   If omitted, the Codex Action selects its default; confirm availability in the
   OpenAI project rather than assuming a hardcoded model will remain available.
5. Start a manual workflow run on the default branch and review the logs/checks
   and resulting draft PR before relying on unattended runs.

A local file is not an active remote schedule. Without the API secret the
workflow exits preflight with a setup message. Without repository authentication
this local development session cannot publish files or configure the secret.

## Schedule and cost

Cron: `17 */2 * * *`, at 00:17, 02:17, …, 22:17 UTC, also the local time in Accra.
GitHub may delay scheduled jobs, and scheduled workflows run from the default
branch only. This is a two-hour maintenance opportunity, not a promise that code
changes every two hours. No justified change means no PR.

API usage is charged to the configured OpenAI project and GitHub runner usage
also applies. Configure project budgets/alerts appropriate to up to 12 daily
opportunities. Each proposal job has a 30-minute limit; verification has a
15-minute limit. Set `CODEX_MAINTENANCE_PAUSED=true` to pause all proposals, or
disable the workflow. Rotate the project key according to your operating policy.

## Execution and privilege separation

1. **Preflight:** read-only GitHub access; require the API secret and check for an
   open `ai/maintenance-*` PR. Skip if one is already awaiting review.
2. **Propose:** check out the exact default-branch SHA without persisted GitHub
   credentials, run the lightweight baseline, then invoke Codex with workspace
   sandboxing and `drop-sudo`. The API key is scoped to the Codex Action step.
   Afterwards only Git diff export and artifact upload occur; no proposed code
   is executed in this job after the agent.
3. **Verify:** use a fresh runner with no OpenAI key or repository write token.
   Run the trusted base version of the patch validator before applying changes,
   then rerun Python/Nextflow/synthetic workflow validation independently.
4. **Publish:** another runner receives GitHub write permissions, rechecks the
   scope and base SHA, applies the verified patch, and creates a unique branch
   and draft PR. It does not execute proposed repository code or merge the PR.

If the default branch changes during validation, publication is deferred to a
new scheduled run. Empty patches do not publish. Any baseline/check/scope failure
blocks publication. GitHub-token-created PRs may not trigger ordinary PR workflows;
the separate verification job supplies the checks before publication. A maintainer
can rerun CI when promoting the draft for review.

## Allowed work and hard limits

The prompt is [maintenance.md](../.github/codex/maintenance.md), governed by
[AGENTS.md](../AGENTS.md). The independent
[patch validator](../scripts/check_maintenance_patch.py) allows only:

- Root README files and Markdown under `docs/`, excluding nested agent policies.
- `bin/validate_inputs.py` and `bin/generate_report.py` for input/presentation fixes.
- Top-level `tests/test_*.py` regression coverage.

It rejects changes to scientific utilities, `.github/`, policies, the validator,
modules/config, references/assets/databases, dependencies and other paths. It
also rejects deletions, renames, binary changes, permission changes, symlinks and
executable additions. A patch may change at most 12 files, 400 added/deleted lines
and 200 KB. These are workflow-enforced limits, not only prompt requests.

The agent must not alter scientific conclusions through report code, weaken
tests or invent results. Scope checks do not replace human review of code meaning.
Scientific tool/database upgrades require a separately reviewed task and full
validation, not a two-hour automatic dependency refresh.

## Troubleshooting

Check preflight for missing key, paused variable or pending PR. Check propose for
API/model permissions and baseline failures. Check verify for rejected paths,
patch conflicts or regressions. Check publish for changed base SHA or disabled PR
creation. Review artifacts within their seven-day retention window. No credentials
or patient/genome datasets belong in those artifacts.
