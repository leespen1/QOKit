# How to Use the Research Log System

## Quick Start

The research log tracks experimental results and open questions for the IEEE
Quantum Week paper. Everything lives in `research_log/`.

### Key Files

| File | Purpose | When to read |
|------|---------|-------------|
| `next_steps.md` | Prioritized research queue (P0/P1/P2/DONE) | Start of any research session |
| `index.md` | One-line summary of every result | Before running an experiment (check novelty) |
| `entries/*.md` | Full write-ups of each experiment | When you need details on a past result |
| `README.md` | Entry template and conventions | When writing a new entry |

### Two Operating Modes

**Active mode** — You're at the computer, guiding the research.

Tell Claude something like:
- "Work on research" — Claude reads next_steps, proposes the top P0 item, waits for approval
- "Run the triangle proxy fitting experiment" — Claude runs a specific item
- "Add a P1 item about testing on regular graphs" — Claude updates next_steps
- "What's in the research log?" — Claude summarizes index.md
- "Reprioritize: move depth scaling to P0" — Claude reshuffles priorities

**Autonomous mode** — Claude works through the queue unattended (e.g., overnight).

Tell Claude:
- "Go autonomous" or "overnight mode"

Claude will:
1. Pick the top P0 item from next_steps.md
2. Run the experiment
3. Log results to entries/, update index.md, update next_steps.md
4. Move to the next P0 item
5. On errors or ambiguity: make a reasonable assumption or skip, then continue
6. Stop when P0 is empty or critically stuck

### Managing the Queue

The priority levels in `next_steps.md` are:

- **P0** — Do next. These are the items Claude picks up in autonomous mode.
- **P1** — Do soon. Promote to P0 when the current P0 batch is done.
- **P2** — Backlog. Ideas worth tracking but not urgent.
- **DONE** — Completed items with links to their log entries.

You can ask Claude to add, remove, reprioritize, or edit items at any time.

### What Gets Logged

Claude logs a result when:
- The experiment produced a significant or novel finding
- The result confirms, refutes, or refines a prior result
- A similar experiment has not already been logged (checked via index.md)

Each entry includes: motivation, setup (script, parameters, output paths),
key findings, significance for the paper, and next steps arising.

### Reviewing Progress

- **Quick scan:** Read `index.md` — one line per result, sortable by date or paper section
- **Deep dive:** Read a specific entry in `entries/`
- **Research direction:** Read `next_steps.md` to see what's done, what's next, and what's queued
