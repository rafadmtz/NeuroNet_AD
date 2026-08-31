#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="$PROJECT_ROOT/output"
TABLA="$PROJECT_ROOT/data/conteocell_todos.tsv"
CSV="$PROJECT_ROOT/data/adstate_fem.csv"

# Lista de subtipos celulares
SUBTIPOS=(
    "Astrocyte" "Chandelier" "Endothelial"
    "L2-3-IT" "L4-IT" "L5-6-NP" "L5-ET" "L5-IT"
    "L6b" "L6-CT" "L6-IT-Car3" "L6-IT"
    "Lamp5-Lhx6" "Lamp5" "Microglia-PVM"
    "Oligodendrocyte" "OPC" "Pax6" "Pvalb"
    "Sncg" "Sst-Chodl" "Sst" "Vip" "VLMC"
)

# Header: ID, AD_State, Sex, y una columna por subtipo
HEADER="ID\tAD_State\tSex"
for subtipo in "${SUBTIPOS[@]}"; do
    HEADER+="\t${subtipo}"
done
printf "%b\n" "$HEADER" > "$TABLA"

# Iterar sobre cada donante en el CSV
sed '1d' "$CSV" | while IFS=',' read -r id ad_state sex || [[ -n "$id" ]]; do
    # Limpiar caracteres raros
    id="${id//$'\r'/}"
    id="${id//\"/}"
    id="$(echo "$id" | xargs)"
    ad_state="${ad_state//$'\r'/}"
    ad_state="${ad_state//\"/}"
    ad_state="$(echo "$ad_state" | xargs)"
    sex="${sex//$'\r'/}"
    sex="${sex//\"/}"
    sex="$(echo "$sex" | xargs)"

    ROW="${id}\t${ad_state}\t${sex}"

    for subtipo in "${SUBTIPOS[@]}"; do
        file=$(find "$DIR" -type f -name "*${id}*_${subtipo}.tsv" | head -n 1)
        if [[ -n "$file" ]]; then
            NUM=$(awk -F'\t' 'NR==1 {print NF-1; exit}' "$file")
        else
            NUM=0
        fi
        ROW+="\t${NUM}"
    done

    printf "%b\n" "$ROW" >> "$TABLA"
done

echo "Listo: $TABLA"
