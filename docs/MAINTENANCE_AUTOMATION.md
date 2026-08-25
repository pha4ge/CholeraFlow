# Weekly OpenAI Codex maintenance

CholeraFlow includes a guarded GitHub Actions workflow at:

```text
.github/workflows/weekly-codex-maintenance.yml
```

It runs weekly and can also be started manually. The workflow uses OpenAI's
official `openai/codex-action@v1` and is intentionally designed to **propose a
pull request, never merge into `main`**.

## One-time GitHub setup

1. Push this repository to GitHub.
2. Open **Settings → Secrets and variables → Actions**.
3. Add a repository secret named `OPENAI_API_KEY` containing an OpenAI API key
   that is allowed to use Codex.
4. In **Settings → Actions → General → Workflow permissions**, allow GitHub
   Actions to create pull requests if your repository/org policy disables it.
5. Keep branch protection/review requirements enabled on `main`.

No API key is stored in the repository. The Codex job has read-only GitHub
content permission and uses the official `:workspace` Codex permission profile
for repository edits. It edits an isolated checkout, runs tests, and exports a
patch. A separate publishing job receives GitHub write permission and creates a
new `ai/weekly-maintenance-*` branch plus pull request. This separation reduces
the privileges available to the coding agent.

## Schedule

The default cron is Monday at 06:00 UTC. Change only the cron expression in the
workflow if a different maintenance window is preferred.

## What it is allowed to do

The persistent policy is in `AGENTS.md` and
`.github/codex/weekly-maintenance.md`. Scientific methodology changes are
explicitly excluded from autonomous maintenance and must be reviewed by a
human.

## Guardrail files

The weekly agent is prevented from editing its own governing files: `AGENTS.md`,
`.github/codex/weekly-maintenance.md`, and all files under `.github/workflows/`.
Changes to those files require a separate human-directed maintenance task.
