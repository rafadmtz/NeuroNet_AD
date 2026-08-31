#!/usr/bin/env bash
# Encadena calculo de MI (entorno con numba_mi) -> analisis de red (entorno
# con igraph), un solo proceso, un solo nohup.
#
# Ajusta ENV_MI y ENV_ANALISIS a los nombres reales de tus entornos conda.
# Se asume que este script vive en la raiz del proyecto y los .py en src/.
#
# Uso:
#   nohup bash run_all.sh > logs/numba_allvip_analisis.log 2>&1 &
set -euo pipefail
 
ENV_MI=numba_mi
ENV_ANALISIS=neuronet
 
source "$(conda info --base)/etc/profile.d/conda.sh"
 
echo "=========== Paso 1/2: calculo de MI ($ENV_MI) ==========="
conda activate "$ENV_MI"
python src/01_calcular_mi.py
conda deactivate
 
echo "=========== Paso 2/2: analisis de red ($ENV_ANALISIS) ==========="
conda activate "$ENV_ANALISIS"
python src/02_analizar_redes.py
conda deactivate
 
echo "Listo."