import matplotlib.pyplot as plt   # Biblioteca de gráficos (plt é o apelido usual)

# ============================================================
# GRÁFICOS DAS CURVAS HIDROSTÁTICAS
# Recebe a lista de resultados de calcular_curva() (um dicionário por calado) e desenha uma
# característica em função do calado. A mesma função serve para qualquer característica:
# basta trocar a "chave" e os textos dos eixos.
# Segue a convenção do IC: o calado fica no eixo vertical.
# ============================================================
#AAAAAAAA
def plotar_curva(curva, chave, titulo, rotulo_x, caminho_saida=None):
    # curva: lista de dicionários, saída de calcular_curva()
    # chave: nome da característica no dicionário (ex.: 'area_molhada')
    # titulo e rotulo_x: textos do título e do eixo horizontal
    # caminho_saida: se informado, também salva o gráfico nesse arquivo (ex.: "area_molhada.png")
    calados = [resultado['calado'] for resultado in curva]   # lista com o calado de cada resultado
    valores = [resultado[chave] for resultado in curva]      # lista com a característica pedida, na mesma ordem

    figura, eixo = plt.subplots(figsize=(7, 6))              # cria a figura (a "folha") e o eixo (o gráfico dentro dela)
    eixo.plot(valores, calados, marker='o', markersize=4)    # x = característica, y = calado; marker mostra os pontos calculados
    eixo.set_xlabel(rotulo_x)                                # texto do eixo horizontal
    eixo.set_ylabel("Calado (m)")                            # texto do eixo vertical
    eixo.set_title(titulo)                                   # título do gráfico
    eixo.grid(True)                                          # grade de fundo, para facilitar a leitura

    if caminho_saida is not None:                            # só salva se um caminho foi dado
        figura.savefig(caminho_saida, dpi=150)

    #plt.show()                                               # abre a janela do gráfico (o programa espera você fechá-la)