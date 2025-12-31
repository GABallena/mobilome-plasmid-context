#!/usr/bin/env bash
# run_mge_oap_local.sh — MGE OAP runner (laptop-friendly, with logs & resume)
# - Uses SHORT-ID MGE DB by default (mge_db.short.fasta)
# - Gentle CPU/IO via nice/ionice
# - Auto-detects fq/fq.gz
# - Resume mode, sample limiting, per-sample logs, failures.txt

set -Eeuo pipefail

# ========= USER TUNABLES (env-overridable) ===================================
READS_DIR="${READS_DIR:-$PWD/trimmed_reads}"            # where *_R1_*trimmed.fq(.gz) live
WORK_BASE="${WORK_BASE:-$PWD/mge_oap_work}"             # outputs here
ARGS_OAP_DB="${ARGS_OAP_DB:-$PWD/MGE_db/mge_db.short.fasta}"  # <<< short FASTA (indexed)

THREADS="${THREADS:-4}"                                  # <= keep small on laptop (2–4)
NICE_LEVEL="${NICE_LEVEL:-12}"                           # 10–15 = gentle
IONICE_CLASS="${IONICE_CLASS:-2}"                        # 2 = best-effort
IONICE_LEVEL="${IONICE_LEVEL:-7}"                        # 0 (high) .. 7 (low)

# Stage filters (can relax/tighten)
FORMAT_EXT="${FORMAT_EXT:-}"                             # auto if empty
E1="${E1:-1e-10}"                                        # stage_one evalue 1
E2="${E2:-3}"                                            # stage_one evalue 2
ID1="${ID1:-45}"                                         # stage_one min identity
QCOV1="${QCOV1:-0}"                                      # stage_one min query coverage

EVAL2="${EVAL2:-1e-7}"                                   # stage_two evalue
PID2="${PID2:-80}"                                       # stage_two identity
QCOV2="${QCOV2:-75}"                                     # stage_two coverage
ALNLEN2="${ALNLEN2:-25}"                                 # stage_two min align length

RESUME_ONLY="${RESUME_ONLY:-0}"                          # 1 = only run missing .stage2.done
MAX_SAMPLES="${MAX_SAMPLES:-0}"                          # 0 = all, else limit count
ULIMIT_NOFILE="${ULIMIT_NOFILE:-2048}"                   # file handles
# ============================================================================

STAGE1_CMD="${STAGE1_CMD:-args_oap stage_one}"
STAGE2_CMD="${STAGE2_CMD:-args_oap stage_two}"

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

# ---------- logging helpers ----------
log_fail_init(){
  mkdir -p "$WORK_BASE"
  FAIL_FILE="$WORK_BASE/failures.txt"
  if [[ ! -f "$FAIL_FILE" ]]; then : > "$FAIL_FILE"; fi
  echo "===== $(ts) run begin =====" >> "$FAIL_FILE"
}
record_failure(){ # $1=sample $2=stage $3=logpath $4=note?
  local s="$1" st="$2" logp="$3" note="${4:-}"
  echo "$(ts)  SAMPLE=${s}  STAGE=${st}  LOG=${logp}  ${note}" >> "$FAIL_FILE"
}

# ---------- helpers ----------
detect_format(){
  if [[ -n "$FORMAT_EXT" ]]; then echo "$FORMAT_EXT"; return; fi
  if ls "$READS_DIR"/*_R1_*trimmed.fq.gz >/dev/null 2>&1; then echo "fq.gz"; return; fi
  if ls "$READS_DIR"/*_R1_*trimmed.fastq.gz >/dev/null 2>&1; then echo "fastq.gz"; return; fi
  if ls "$READS_DIR"/*_R1_*trimmed.fq >/dev/null 2>&1; then echo "fq"; return; fi
  echo "fastq"
}
pair_for_r1(){ # prints R2 or fails
  local r1="$1" r2=""
  r2="${r1/_R1_/_R2_}"; [[ -f "$r2" ]] || r2="${r1/_R1/_R2}"
  [[ -f "$r2" ]] || return 1
  printf '%s\n' "$r2"
}
have_prog(){ command -v "$1" >/dev/null 2>&1; }

ensure_blast_db(){
  # Build BLAST indexes if missing (nin/nsq/nhr). Uses makeblastdb (nucl).
  local f base
  f="${1-}"
  if [[ -z "$f" ]]; then
    echo "[WARN] ensure_blast_db: missing fasta path; skipping index check/build."
    return 0
  fi
  base="${f%.*}"
  # Detect either single-file indexes (base.nin/nsq/nhr) or split volumes (base.00.nhr etc.)
  if compgen -G "${base}.nin" >/dev/null ||
     compgen -G "${base}.nsq" >/dev/null ||
     compgen -G "${base}.nhr" >/dev/null ||
     compgen -G "${base}."*".nin" >/dev/null ||
     compgen -G "${base}."*".nsq" >/dev/null ||
     compgen -G "${base}."*".nhr" >/dev/null ||
     compgen -G "${base}.nal" >/dev/null ||
     compgen -G "${base}.ndb" >/dev/null ||
     compgen -G "${base}.not" >/dev/null ||
     compgen -G "${base}.ntf" >/dev/null ||
     compgen -G "${base}.nto" >/dev/null; then
    return 0
  fi
  if ! have_prog makeblastdb; then
    echo "[WARN] makeblastdb not found; continuing (assuming prebuilt indexes)."
    return 0
  fi
  echo "[INFO] Building BLAST DB indexes for $f ..."
  makeblastdb -in "$f" -dbtype nucl -out "$base" >/dev/null 2>&1 || {
    echo "[ERROR] makeblastdb failed for $f"; return 1; }
}

# ---------- main steps ----------
preflight(){
  [[ -d "$READS_DIR" ]] || { echo "ERROR: Reads dir not found: $READS_DIR" >&2; exit 1; }
  [[ -f "$ARGS_OAP_DB" ]] || { echo "ERROR: DB fasta not found: $ARGS_OAP_DB" >&2; exit 1; }
  mkdir -p "$WORK_BASE/sample_in" "$WORK_BASE/samples" "$WORK_BASE/tmp"
  : "${TMPDIR:=$WORK_BASE/tmp}"; export TMPDIR
  ulimit -n "$ULIMIT_NOFILE" || true
  ensure_blast_db "$ARGS_OAP_DB" || echo "[WARN] Could not ensure BLAST DB — continuing."
}

collect_r1s(){
  mapfile -t R1S < <(find "$READS_DIR" -maxdepth 1 -type f -name "*_R1_*trimmed.fq.gz" | sort)
  if [[ ${#R1S[@]} -eq 0 ]]; then
    mapfile -t R1S < <(find "$READS_DIR" -maxdepth 1 -type f -name "*_R1_*trimmed.fq" | sort)
  fi
  echo "Found ${#R1S[@]} R1 files."
  if [[ "$RESUME_ONLY" == "1" ]]; then
    local tmp=()
    for r1 in "${R1S[@]}"; do
      local b; b="$(basename "$r1")"; b="${b%%_R1_*}"
      [[ -f "$WORK_BASE/samples/$b/.stage2.done" ]] || tmp+=("$r1")
    done
    R1S=("${tmp[@]}")
    echo "Resuming: ${#R1S[@]} to run."
  fi
  if (( MAX_SAMPLES>0 )) && (( ${#R1S[@]}>MAX_SAMPLES )); then
    R1S=("${R1S[@]:0:MAX_SAMPLES}")
    echo "Limiting to first $MAX_SAMPLES samples."
  fi
}

run_one(){
  local r1="$1" base r2 s_in s_out log done1 done2 fmt rc
  base="$(basename "$r1")"; base="${base%%_R1_*}"
  r2="$(pair_for_r1 "$r1")" || { echo "[FAIL] $base pair"; record_failure "$base" "pair" "(none)" "R2 not found"; return 1; }

  s_in="${WORK_BASE}/sample_in/${base}"
  s_out="${WORK_BASE}/samples/${base}"
  log="${s_out}/logs/${base}.log"
  done1="${s_out}/.stage1.done"; done2="${s_out}/.stage2.done"
  mkdir -p "$s_in" "$s_out/logs"

  # Symlink inputs (keeps your original pattern)
  ln -sf "$r1" "$s_in/"; ln -sf "$r2" "$s_in/"

  fmt="$(detect_format)"

  echo "== [$base] =="

  # ---- Stage 1 ----
  if [[ -f "$done1" && -s "${s_out}/metadata.txt" && -s "${s_out}/extracted.fa" ]]; then
    echo "[SKIP] Stage 1"
  else
    echo "[RUN ] Stage 1..."
    {
      echo "== $(ts) :: Stage 1 =="
      echo "DB: $ARGS_OAP_DB"
      echo "FMT: $fmt  THREADS: $THREADS  nice: $NICE_LEVEL  ionice: $IONICE_CLASS/$IONICE_LEVEL"
      echo "CMD: $STAGE1_CMD -i $s_in -o $s_out -t $THREADS -f $fmt --e1 $E1 --e2 $E2 --id $ID1 --qcov $QCOV1"
    } >> "$log"
    set +e
    OMP_NUM_THREADS="$THREADS" BLAST_USAGE_THREADS="$THREADS" ARGS_OAP_DB="$ARGS_OAP_DB" \
    nice -n "$NICE_LEVEL" ionice -c "$IONICE_CLASS" -n "$IONICE_LEVEL" \
    $STAGE1_CMD -i "$s_in" -o "$s_out" -t "$THREADS" -f "$fmt" \
      --e1 "$E1" --e2 "$E2" --id "$ID1" --qcov "$QCOV1" >> "$log" 2>&1
    rc=$?; set -e
    if [[ $rc -ne 0 || ! -s "${s_out}/metadata.txt" || ! -s "${s_out}/extracted.fa" ]]; then
      echo "[FAIL] $base (Stage 1). See $log"; record_failure "$base" "stage1" "$log"; return 1
    fi
    touch "$done1"
  fi

  # ---- Stage 2 ----
  if [[ -f "$done2" ]]; then
    echo "[SKIP] Stage 2"
    echo "[OK  ] $base"
    return 0
  else
    echo "[RUN ] Stage 2..."
    {
      echo "== $(ts) :: Stage 2 =="
      echo "CMD: $STAGE2_CMD -i $s_out -t $THREADS --e $EVAL2 --id $PID2 --qcov $QCOV2 --length $ALNLEN2"
    } >> "$log"
    set +e
    OMP_NUM_THREADS="$THREADS" BLAST_USAGE_THREADS="$THREADS" ARGS_OAP_DB="$ARGS_OAP_DB" \
    nice -n "$NICE_LEVEL" ionice -c "$IONICE_CLASS" -n "$IONICE_LEVEL" \
    $STAGE2_CMD -i "$s_out" -t "$THREADS" \
      --e "$EVAL2" --id "$PID2" --qcov "$QCOV2" --length "$ALNLEN2" >> "$log" 2>&1
    rc=$?; set -e
    if [[ $rc -ne 0 ]]; then
      echo "[FAIL] $base (Stage 2). See $log"; record_failure "$base" "stage2" "$log"; return 1
    fi
    touch "$done2"
  fi

  echo "[OK  ] $base"
}

main(){
  preflight
  log_fail_init

  collect_r1s
  echo "Running with THREADS=$THREADS  (nice=$NICE_LEVEL, ionice=$IONICE_CLASS/$IONICE_LEVEL)"
  for r1 in "${R1S[@]}"; do run_one "$r1" || true; done

  echo "----- SUMMARY -----"
  echo "Inputs    : ${#R1S[@]}"
  echo "Outputs in: $WORK_BASE/samples/"
  echo "Failures  : see $WORK_BASE/failures.txt (if any)"
}
main "$@"
