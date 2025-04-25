from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
import numpy
import copy
from fastaReader import fastaReader

def run_bfoa(
    numeroDeBacterias=5,
    numRandomBacteria=2,
    iteraciones=5,
    tumbo=2,
    nado=3,
    dAttr=0.1,
    wAttr=0.002,
    hRep=0.1,
    wRep=0.001,
    max_processes=4
):
    """
    Ejecuta el algoritmo parallel BFOA y devuelve:
      - veryBestFitness: mejor fitness obtenido
      - globalNFE: total de evaluaciones de función objetivo
      - veryBestBlosum: puntuación BLOSUM del mejor individuo

    Parámetros no activos:
      - numRandomBacteria: número de bacterias iniciales aleatorias (no implementado internamente)
    """
    # Cargar secuencias
    secuencias = fastaReader().seqs
    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])

    # Inicializaciones
    globalNFE = 0
    manager = Manager()
    numSec = len(secuencias)
    poblacion = manager.list(range(numeroDeBacterias))

    # Población inicial: copia directa de las secuencias
    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)

    operadorBacterial = bacteria(numeroDeBacterias, max_processes=max_processes)
    veryBestIdx = None
    veryBestFitness = None
    veryBestBlosum = None

    start_time = time.time()
    poblacionInicial()

    # Evolución
    for it in range(iteraciones):
        # 1) Tumble + registro de movimientos
        moves = operadorBacterial.tumbo(numSec, poblacion, tumbo)

        # 2) Evaluación tras tumble
        operadorBacterial.cuadra(numSec, poblacion)
        operadorBacterial.creaGranListaPares(poblacion)
        operadorBacterial.evaluaBlosum()
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        operadorBacterial.creaTablaInteraction()
        operadorBacterial.creaTablaFitness()

        # Acumular NFE tras evaluación
        globalNFE += operadorBacterial.getNFE()

        best_local = operadorBacterial.obtieneBest(globalNFE)[1]

        # 3) Swim: reaplicar moves hasta nado pasos mientras mejore
        pasos_swim = nado
        while pasos_swim > 0:
            operadorBacterial.tumbo_apply_moves(numSec, poblacion, moves)
            operadorBacterial.cuadra(numSec, poblacion)
            operadorBacterial.creaGranListaPares(poblacion)
            operadorBacterial.evaluaBlosum()
            operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
            operadorBacterial.creaTablaInteraction()
            operadorBacterial.creaTablaFitness()

            # Acumular NFE de swim
            globalNFE += operadorBacterial.getNFE()

            nuevo_best = operadorBacterial.obtieneBest(globalNFE)[1]

            if nuevo_best > best_local:
                best_local = nuevo_best
                pasos_swim -= 1
            else:
                break

        # 4) Actualizar global best y reemplazar peor
        if veryBestFitness is None or best_local > veryBestFitness:
            veryBestFitness = best_local
            veryBestIdx     = operadorBacterial.obtieneBest(globalNFE)[0]
            veryBestBlosum  = operadorBacterial.blosumScore[veryBestIdx]

        operadorBacterial.replaceWorst(poblacion, veryBestIdx)
        operadorBacterial.resetListas(numeroDeBacterias)

    # Retornar métricas
    return veryBestFitness, globalNFE, veryBestBlosum


if __name__ == "__main__":
    # Ejemplo de uso: 3 iteraciones internas y swim de 3 pasos
    fitness, interacciones, blosum = run_bfoa(
        numeroDeBacterias=4,
        numRandomBacteria=1,
        iteraciones=3,
        tumbo=2,
        nado=3,
        max_processes=2
    )
    print(f"\nVery Best Fitness: {fitness}")
    print(f"Interacciones: {interacciones}")
    print(f"BLOSUM: {blosum}")
