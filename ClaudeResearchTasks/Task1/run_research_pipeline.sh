#!/usr/bin/env bash
set -euo pipefail

TASK_FILE="research_task1.md"
DRAFT_FILE="research_task1_draft.md"
CRITIQUE_FILE="research_task1_critique.md"
FINAL_FILE="research_task1_final.md"
LOG_FILE="research_task1_pipeline.log"

{
  echo "=== Stage 1: Draft ==="
  claude --dangerously-skip-permissions -p < "$TASK_FILE" > "$DRAFT_FILE"

  echo "=== Stage 2: Critique ==="
  cat <<EOF | claude --dangerously-skip-permissions -p > "$CRITIQUE_FILE"
You are reviewing a draft produced for the following task.

First, read the original task:
----- TASK -----
$(cat "$TASK_FILE")
----- END TASK -----

Now critique this draft harshly and concretely.
Focus on:
1. missing requirements
2. weak reasoning
3. unsupported claims
4. vagueness
5. bad structure
6. places where the answer likely stopped too early

Return:
- a short overall verdict
- a numbered list of problems
- a numbered list of exact revision actions

----- DRAFT -----
$(cat "$DRAFT_FILE")
----- END DRAFT -----
EOF

  echo "=== Stage 3: Revise ==="
  cat <<EOF | claude --dangerously-skip-permissions -p > "$FINAL_FILE"
Revise the draft below using the critique.

Requirements:
- satisfy the original task fully
- fix every critique item
- do not merely comment on the draft; rewrite it into a stronger final deliverable
- make reasonable assumptions instead of asking questions
- output only the final revised result

----- ORIGINAL TASK -----
$(cat "$TASK_FILE")
----- END ORIGINAL TASK -----

----- CRITIQUE -----
$(cat "$CRITIQUE_FILE")
----- END CRITIQUE -----

----- DRAFT -----
$(cat "$DRAFT_FILE")
----- END DRAFT -----
EOF

  echo "Pipeline complete."
} 2>&1 | tee "$LOG_FILE"
