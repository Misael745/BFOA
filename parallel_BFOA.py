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
      - nado: longitud de "nado" tras cada tumble (no implementado internamente)
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
        # fase tumble
        operadorBacterial.tumbo(numSec, poblacion, tumbo)
        # (no hay fase de nado implementada)
        operadorBacterial.cuadra(numSec, poblacion)
        operadorBacterial.creaGranListaPares(poblacion)
        operadorBacterial.evaluaBlosum()
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        operadorBacterial.creaTablaInteraction()
        operadorBacterial.creaTablaFitness()

        globalNFE += operadorBacterial.getNFE()
        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
        currentBlosum = operadorBacterial.blosumScore[bestIdx]

        # Actualizar global best
        if veryBestFitness is None or bestFitness > veryBestFitness:
            veryBestIdx = bestIdx
            veryBestFitness = bestFitness
            veryBestBlosum = currentBlosum

        operadorBacterial.replaceWorst(poblacion, veryBestIdx)
        operadorBacterial.resetListas(numeroDeBacterias)

    # Retornar métricas
    return veryBestFitness, globalNFE, veryBestBlosum

if __name__ == "__main__":
    # Ejemplo de uso: 3 iteraciones internas, 1 bacteria aleatoria (no activa)
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
