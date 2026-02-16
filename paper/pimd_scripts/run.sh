#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 # crashes sometimes if I dont have this
LMP="/home/emeitz/software/lammps_ipi/build/lmp"
LMP_TEMPLATE="in.ipi_lj.lmp"
IPI_TEMPLATE="neon_cv.xml"


# -------------------
# Temperatures (K) as integers
# -------------------
TEMPS_K=(4 6 8 10 12 14 16 18 20 22 24)

# -------------------
# i-PI / PIMD controls
# -------------------
DT_FS="3.0"
TAU_FS="100.0"
DATA_STRIDE="10"
MAX_CORES=40

# -------------------
# LJ parameters (Ne)
# -------------------
EPS="0.0032135"   # eV
SIG="2.782"       # Å
RC=$(python - <<PY
sig=float("$SIG")
print(2.5*sig)
PY
)


gen_seed () {
  python - <<'PY'
import secrets
print(secrets.randbelow(2**31 - 1) + 1)
PY
}

beads_for_T () {
  python - <<PY
import math
T=float("$1")
print(int(math.ceil(550.0/T)))
PY
}

# -------------------
# Compute NBEADS for all temperatures
# -------------------
NBEADS_LIST=()
for T in "${TEMPS_K[@]}"; do
  NBEADS_LIST+=("$(beads_for_T "$T")")
done

# -------------------
# Allocate NCLIENT per temperature (bead-weighted, integer, capped)
# Uses a "largest remainder" apportionment + cap redistribution.
# -------------------
NCLIENT_LIST=($(python3 - "$MAX_CORES" "${TEMPS_K[@]}" -- "${NBEADS_LIST[@]}" <<'PY'
import sys, math

MAX_CORES = int(sys.argv[1])

# argv layout: MAX_CORES, temps..., --, beads...
args = sys.argv[2:]
sep = args.index("--")
temps = list(map(int, args[:sep]))
beads = list(map(int, args[sep+1:]))

n = len(beads)
assert n > 0
assert len(temps) == n

# At least 1 client per sim (cap by beads)
clients = [min(1, b) for b in beads]

remaining = MAX_CORES - sum(clients)
if remaining < 0:
    remaining = 0

while remaining > 0:
    room = [b - c for b, c in zip(beads, clients)]
    active = [i for i, r in enumerate(room) if r > 0]
    if not active:
        break

    weights = [beads[i] for i in active]
    wsum = sum(weights)
    if wsum <= 0:
        break

    frac = [remaining * (w / wsum) for w in weights]
    extra = [int(math.floor(x)) for x in frac]
    used = sum(extra)

    if used == 0:
        i = max(active, key=lambda j: beads[j])
        clients[i] += 1
        remaining -= 1
        continue

    for idx, e in zip(active, extra):
        add = min(e, beads[idx] - clients[idx])
        clients[idx] += add

    remaining = MAX_CORES - sum(clients)

    if remaining <= 0:
        break

    room = [b - c for b, c in zip(beads, clients)]
    active = [i for i, r in enumerate(room) if r > 0]
    if not active:
        break

    weights = [beads[i] for i in active]
    wsum = sum(weights)
    if wsum <= 0:
        break

    frac = [(MAX_CORES - sum(clients)) * (w / wsum) for w in weights]
    order = sorted(
        range(len(active)),
        key=lambda k: frac[k] - math.floor(frac[k]),
        reverse=True
    )

    for k in order:
        if remaining <= 0:
            break
        i = active[k]
        if clients[i] < beads[i]:
            clients[i] += 1
            remaining -= 1

    remaining = MAX_CORES - sum(clients)

clients = [min(c, b) for c, b in zip(clients, beads)]

for c in clients:
    print(c)
PY
))

SUMCLIENTS=0
for x in "${NCLIENT_LIST[@]}"; do
  SUMCLIENTS=$((SUMCLIENTS + x))
done

if (( SUMCLIENTS != MAX_CORES )); then
  echo "WARNING: SumClients=${SUMCLIENTS} != MAX_CORES=${MAX_CORES}"
  echo "This can happen if all sims hit the cap NCLIENTS=NBEADS."
fi

echo "===================================================="
echo "MAX_CORES = $MAX_CORES"
echo "Temps     = ${TEMPS_K[*]}"
echo "Beads     = ${NBEADS_LIST[*]}"
echo "Clients   = ${NCLIENT_LIST[*]}"
echo "SumClients= ${SUMCLIENTS}"
echo "===================================================="

# -------------------
# Run ALL temperatures concurrently (each has its own unix socket address)
# -------------------
pids_all=()

for i in "${!TEMPS_K[@]}"; do
  TEMP_K="${TEMPS_K[$i]}"
  NBEADS="${NBEADS_LIST[$i]}"
  NCLIENTS="${NCLIENT_LIST[$i]}"
  seed="$(gen_seed)"

  d="run_T${TEMP_K}K"
  mkdir -p "$d"

  # To initialize LAMMPS 
  DATA_NAME="neon_T${TEMP_K}.data"
  DATA_FILE="./structures/lammps_data/${DATA_NAME}"
  if [[ ! -f "$DATA_FILE" ]]; then
    echo "ERROR: missing structure file: $DATA_FILE"
    exit 1
  fi

  # To initialize i-pi
  XYZ_NAME="neon_T${TEMP_K}.xyz"
  XYZ_FILE="./structures/xyz/${XYZ_NAME}"
  if [[ ! -f "$XYZ_FILE" ]]; then
    echo "ERROR: missing structure file: $XYZ_FILE"
    exit 1
  fi

  cp "$DATA_FILE" "$d"/
  cp "$XYZ_FILE" "$d"/
  cp "$LMP_TEMPLATE" "$d"/in.lmp
  cp "$IPI_TEMPLATE" "$d"/input.xml

  # unique socket name per temperature run
  ADDR="lj_T${TEMP_K}K"
  SOCKET="/tmp/ipi_${ADDR}"

  sed -i \
    -e "s/NBEADS/${NBEADS}/g" \
    -e "s/NCLIENTS/${NCLIENTS}/g" \
    -e "s/TEMP_KK/${TEMP_K}K/g" \
    -e "s/TEMP_K/${TEMP_K}/g" \
    -e "s/DT_FS/${DT_FS}/g" \
    -e "s/TAU_FS/${TAU_FS}/g" \
    -e "s/DATA_STRIDE/${DATA_STRIDE}/g" \
    -e "s/RNG_SEED/${seed}/g" \
    -e "s/ADDR/${ADDR}/g" \
    "$d"/input.xml

  (
    cd "$d"

    echo "=============================================="
    echo "Running T=${TEMP_K} K | beads=${NBEADS} | clients=${NCLIENTS} | seed=${seed}"
    echo "socket: ${SOCKET}"
    echo "=============================================="

    rm -f "${SOCKET}"

    # start i-PI
    PYTHONUNBUFFERED=1 i-pi input.xml > ipi.log 2>&1 &
    ipipid=$!

    # wait for this temperature's unix socket
    while [[ ! -S "${SOCKET}" ]]; do
      if ! kill -0 "$ipipid" 2>/dev/null; then
        echo "ERROR: i-PI died before creating ${SOCKET}. Check ipi.log"
        exit 1
      fi
      sleep 0.05
    done

    # start LAMMPS clients for THIS temperature
    pids=()
    for ((c=0; c<NCLIENTS; c++)); do
        ${LMP} \
        -in "in.lmp" \
        -screen "none" \
        -var EPS "$EPS" \
        -var SIG "$SIG" \
        -var RC  "$RC" \
        -var ADDR "$ADDR" \
        -var DATA_FILE "$DATA_NAME"
        > "lmp_${c}.log" 2>&1 &
      pids+=($!)
    done

    # cleanup if something crashes
    trap 'kill "$ipipid" "${pids[@]}" 2>/dev/null || true' EXIT

    wait "$ipipid" || true
    for pid in "${pids[@]}"; do
      wait "$pid" || true
    done

    echo "Done T=${TEMP_K} K"
  ) &

  pids_all+=($!)
done

# Wait for all temps
for pid in "${pids_all[@]}"; do
  wait "$pid" || true
done

echo "All temperatures complete."
