import time
import csv
from parallel_BFOA import run_bfoa


def experimento_bfoa(
    n_corridas=30,
    tumbo=2,
    nado=5,
    iteraciones=10,
    archivo_salida="resultados.csv"
):


    with open(archivo_salida, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Corrida", "Fitness", "Tiempo (s)", "Interacciones", "BLOSUM"])

        for i in range(1, n_corridas + 1):
            print(f"\n✨ Ejecutando corrida {i} de {n_corridas}...")
            start = time.time()
            fitness, interacciones, blosum = run_bfoa(
                iteraciones=iteraciones,
                tumbo=tumbo,
                nado=nado
            )
            tiempo = time.time() - start
            writer.writerow([i, fitness, f"{tiempo:.2f}", interacciones, blosum])
            print(f"✅ Corrida {i}: Fitness={fitness}, Tiempo={tiempo:.2f}s, Interacciones={interacciones}, BLOSUM={blosum}")

if __name__ == "__main__":
    experimento_bfoa()
