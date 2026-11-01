from leitura import ler_tabela_cotas  
from splines import (construir_linha_baliza, construir_linha_dagua,
                        montar_malha_final, formar_paineis, vetor_normal_painel, 
                        gerar_grid_com_n_pontos, espelhar_paineis, fechar_fundo,
                        fechar_pontas, indice_coluna_valida)  
from visualizar import gerar_visualizacao_3d
from graficos import plotar_curva
from hidrostatica import paineis_para_array, calcular_hidrostatica, calcular_curva
import numpy as np  
import matplotlib.pyplot as plt                                         
#A parte acima importa bibliotecas do python e funções dos outros código/arquivos dessa pasta

# ============================================================
# Leitura do arquivo do casco selecionado.
# Cada casco tem UM arquivo CSV com balizas, linhas d'agua e cotas juntas (tabela de cotas).
# É chamada a função "ler_tabela_cotas()", e os dados são enviados para a variável "tabela_cotas".
# Só uma embarcação deve estar ativa por vez: as demais ficam comentadas (#).
# Para trocar de casco, comente o bloco ativo e descomente o desejado.
# ============================================================

#Embarcação série 60. Coeficiente de bloco: 0.60
# Argumento: caminho do arquivo do casco.
# O prefixo r"..." indica "raw string", para o Python não interpretar as barras invertidas do caminho do Windows.
tabela_cotas = ler_tabela_cotas(r"C:\Users\Vieira\Documents\Casco-060-4210\Casco-4210.csv")

#Embarcação série 60. Coeficiente de bloco: 0.65
#tabela_cotas = ler_tabela_cotas(r"C:\Users\Vieira\Documents\Casco-065-4211\Casco-4211.csv")

#Embarcação série 60. Coeficiente de bloco: 0.70
#tabela_cotas = ler_tabela_cotas(r"C:\Users\Vieira\Documents\Casco-070-4212\Casco-4212.csv")

#Embarcação série 60. Coeficiente de bloco: 0.75
#tabela_cotas = ler_tabela_cotas(r"C:\Users\Vieira\Documents\Casco-075-4213\Casco-4213.csv")

#Embarcação do livro PNA
#tabela_cotas = ler_tabela_cotas(r"C:\Users\Vieira\Documents\Casco-PNA.csv")

# ============================================================
# Resolução da malha.
# Escolha aqui quantas balizas (posições x) e quantas linhas d'água (alturas z) a malha final terá.
# Os números valem para a malha FINAL: pedir 50 balizas gera exatamente 50 balizas.
# Eles não podem ser menores que a quantidade de pontos reais do casco (balizas com dado e
# alturas z distintas), e todos os pontos reais são mantidos na malha.
# Mais balizas e linhas d'água = mais painéis (número de painéis ~ n_balizas x n_linhas_dagua).
# ============================================================

n_balizas = 73          # número de balizas (posições x) da malha final
n_linhas_dagua = 71    # número de linhas d'água (alturas z) da malha final

# ============================================================
# Grid vertical (alturas z).
# Cada elemento de "tabela_cotas" é uma baliza: baliza[0] é a posição x e baliza[1] é a lista de alturas z.
# Aqui juntamos todas as alturas reais e criamos um grid com exatamente "n_linhas_dagua" alturas,
# que será usado para avaliar as splines das balizas.
# ============================================================

# Lista única com todas as alturas z medidas em todas as balizas
zs_reais = [z for baliza in tabela_cotas for z in baliza[1]]
print(f"Linhas d'água: {n_linhas_dagua} pedidas (mínimo para este casco: {len(set(zs_reais))}).")

# Distribui as linhas d'água pedidas entre as alturas reais, que são todas mantidas
zs_desejados = gerar_grid_com_n_pontos(zs_reais, n_linhas_dagua)

# ============================================================
# Construção das splines das balizas.
# Para cada baliza (seção transversal), ajusta-se uma spline y(z) e avalia-se
# nas alturas refinadas, obtendo uma curva suave por seção.
# ============================================================

# Uma linha de baliza (curva suave) para cada baliza da tabela, avaliada em zs_desejados
malha_balizas = [construir_linha_baliza(baliza, zs_desejados) for baliza in tabela_cotas]

# ============================================================
# Grid longitudinal (posições x).
# Define em quais posições ao longo do comprimento do casco a malha final será calculada.
# Assim como no grid vertical, o grid tem exatamente "n_balizas" posições e mantém as posições x
# das balizas reais (nenhum ponto real é perdido); os pontos extras são distribuídos entre
# as balizas reais de forma que o espaçamento fique o mais uniforme possível.
# ============================================================

xs_balizas = [baliza[0] for baliza in tabela_cotas]  # posição x de cada baliza original
print(f"Balizas: {n_balizas} pedidas (mínimo para este casco: {len(set(xs_balizas))}).")

# Distribui as balizas pedidas entre as balizas reais, que são todas mantidas
xs_desejados = gerar_grid_com_n_pontos(xs_balizas, n_balizas)

# ============================================================
# Construção das splines das linhas d'água e montagem da malha.
# Para cada altura z refinada (índice j), ajusta-se uma spline ao longo de x usando as
# balizas, e avalia-se nas posições xs_desejados.
# Depois, tudo é organizado em uma malha final de pontos 3D.
# ============================================================

# Uma linha d'água (curva ao longo de x) para cada altura z refinada
malha_linhas_dagua = [construir_linha_dagua(malha_balizas, j, xs_desejados)
                       for j in range(len(zs_desejados))]

# Organiza as linhas d'água em uma malha final (pontos organizados por baliza x altura)
malha_final = montar_malha_final(malha_linhas_dagua, len(xs_desejados))

# ============================================================
# Geração dos painéis do casco.
# Os pontos da malha são agrupados em painéis (quadriláteros) que formam a superfície.
# Gera-se um lado (boreste), espelha-se para o outro (bombordo) e fecham-se
# o fundo e as pontas (proa e popa).
# ============================================================

paineis = formar_paineis(malha_final)                 # painéis do lado de boreste
paineis_bombordo = espelhar_paineis(paineis)          # painéis de bombordo, espelhados de boreste
paineis_fundo = fechar_fundo(malha_final)             # painéis que fecham o fundo (quilha)
# Painéis que fecham proa e popa; também retorna os índices das balizas onde o fechamento foi feito
paineis_pontas, idx_proa, idx_popa = fechar_pontas(malha_final)

# Se o fechamento da popa não ocorreu na última baliza, os dados não cobrem todo o comprimento.
# Avisa o usuário e informa quanto do casco ficou sem representação.
if idx_popa != len(malha_final) - 1:
    print(f"Aviso: dado insuficiente na popa. Fechamento feito na baliza x={xs_desejados[idx_popa]:.2f} m, "
          f"em vez de x={xs_desejados[-1]:.2f} m (faltam {xs_desejados[-1]-xs_desejados[idx_popa]:.2f} m sem representar).")

# Mesma verificação para a proa: se o fechamento não ocorreu na primeira baliza, há falta de dados
if idx_proa != 0:
    print(f"Aviso: dado insuficiente na proa. Fechamento feito na baliza x={xs_desejados[idx_proa]:.2f} m.")

# ============================================================
# Lista completa de painéis.
# Junta todos os painéis em uma única lista, calcula o vetor normal de cada um
# e gera a visualização 3D interativa em HTML.
# ============================================================

# Lista completa de painéis: boreste + bombordo + fundo + pontas
todos_paineis = paineis + paineis_bombordo + paineis_fundo + paineis_pontas

# ============================================================
# Cálculos hidrostáticos.
# Os painéis viram um array numpy (o "original", que nunca é alterado) e a meia-nau é definida.
# Cada calado é preparado a partir desse original, em hidrostatica.py.
# Calado específico = lista com um só calado; curva = lista com vários calados (sem if/else).
# Obs.: o maior z da tabela do casco 4210 é o pontal (32 ft), não o calado de projeto (21,33 ft).
# ============================================================

paineis_array = paineis_para_array(todos_paineis)        # painéis originais, formato (N, 4, 3)
x_meia_nau = (xs_balizas[0] + xs_balizas[-1]) / 2        # meia-nau = metade do comprimento coberto pela tabela (x = 0 na proa)

# --- Calado específico ---
calado_unico = 6.5024                                    # calado em metros (aqui, o calado de projeto da Série 60: B/T = 2,5)
resultado = calcular_hidrostatica(paineis_array, calado_unico, x_meia_nau)   # todas as características para este calado
print(f"Calado {resultado['calado']:.3f} m: volume = {resultado['volume']:.2f} m3, "
      f"área molhada = {resultado['area_molhada']:.2f} m2, "
      f"A_WL = {resultado['area_plano_flutuacao']:.2f} m2")

# Posição (índice) das balizas e das linhas d'água REAIS dentro da malha final.
# O visualizador usa isso para desenhar, por padrão, só as linhas reais da tabela de cotas
# (as demais, geradas pelas splines, aparecem com a opção "Todas as linhas").
indices_balizas_reais = [xs_desejados.index(x) for x in xs_balizas]
indices_linhas_reais = [zs_desejados.index(z) for z in sorted(set(zs_reais))]

# --- Curvas hidrostáticas (descomente para gerar os gráficos) ---
calados = np.concatenate([[0.01, 0.1, 0.25, 0.4, 0.6, 0.8], np.linspace(1.0, zs_desejados[-1], 30)])   # 6 calados extras perto da quilha + 30 uniformes de 1 m até o pontal
curva = calcular_curva(paineis_array, calados, x_meia_nau)             # um resultado por calado
plotar_curva(curva, 'area_molhada', "Área molhada x calado (casco 4210)", "Área molhada (m²)")
plotar_curva(curva, 'volume', "Volume x calado (casco 4210)", "Volume (m³)")
plotar_curva(curva, 'area_plano_flutuacao', "Área do plano de flutuação x calado (casco 4210)", "A_WL (m²)")

plt.show()                                                             # abre todos os gráficos de uma vez


# Gera o arquivo HTML com o casco em 3D; os painéis de fechamento (fundo e pontas) são agrupados
#gerar_visualizacao_3d(
#    malha_final,                                          # malha de pontos do casco
#    paineis_boreste=paineis,                              # painéis de boreste
#    paineis_bombordo=paineis_bombordo,                    # painéis de bombordo
#    paineis_fechamento=paineis_fundo + paineis_pontas,    # painéis de fundo, proa e popa
#    caminho_saida="visualizacao_casco.html",              # arquivo de saída
#    indices_balizas_reais=indices_balizas_reais,          # índices das balizas reais na malha
#    indices_linhas_dagua_reais=indices_linhas_reais,      # índices das linhas d'água reais na malha
#)
#AAAAAAAAAA