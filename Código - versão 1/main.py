import numpy as np  
from leitura import ler_tabela_cotas  
from splines import (construir_linha_baliza, construir_linha_dagua,
                        montar_malha_final, formar_paineis, vetor_normal_painel, 
                        gerar_grid_refinado, espelhar_paineis, fechar_fundo,
                        fechar_pontas, indice_coluna_valida)  
from visualizar import gerar_visualizacao_3d  
#A parte acima importa bibliotecas do python e funções dos outros código/arquivos dessa pasta

# ============================================================
# Leitura de arquivos, balizas, linhas d'agua e cotas das embarcações selecionadas.
# É chamado a função "ler_tabela_cotas()", os dados são enviados para a variavel "tabela_cotas".
# Só uma embarcação deve estar ativa por vez: as demais ficam comentadas (#).
# Para trocar de casco, comente o bloco ativo e descomente o desejado.
# ============================================================

#Embarcação série 60. Coeficiente de bloco: 0.60
# Argumentos: (arquivo das balizas, arquivo das linhas d'água, arquivo das cotas).
# O prefixo r"..." indica "raw string", para o Python não interpretar as barras invertidas do caminho do Windows.
tabela_cotas = ler_tabela_cotas(
    r"C:\Users\Vieira\Documents\Casco-060-4210\Balizas-4210.csv",       
    r"C:\Users\Vieira\Documents\Casco-060-4210\Linhas_dagua-4210.csv",  
    r"C:\Users\Vieira\Documents\Casco-060-4210\Cotas-4210.csv",      
)

#Embarcação série 60. Coeficiente de bloco: 0.65
#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-065-4211\Balizas-4211.csv",
#    r"C:\Users\Vieira\Documents\Casco-065-4211\Linhas_dagua-4211.csv",
#    r"C:\Users\Vieira\Documents\Casco-065-4211\Cotas-4211.csv",
#)

#Embarcação série 60. Coeficiente de bloco: 0.70
#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-070-4212\Balizas-4212.csv",
#    r"C:\Users\Vieira\Documents\Casco-070-4212\Linhas_dagua-4212.csv",
#    r"C:\Users\Vieira\Documents\Casco-070-4212\Cotas-4212.csv",
#)

#Embarcação série 60. Coeficiente de bloco: 0.75
#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-075-4213\Balizas-4213.csv",
#    r"C:\Users\Vieira\Documents\Casco-075-4213\Linhas_dagua-4213.csv",
#    r"C:\Users\Vieira\Documents\Casco-075-4213\Cotas-4213.csv",
#)

#Embarcação do livro PNA
#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\balizas.csv",
#    r"C:\Users\Vieira\Documents\linhas_dagua.csv",
#    r"C:\Users\Vieira\Documents\cotas.csv"
#)

# ============================================================
# Refinamento do grid vertical (alturas z).
# Cada elemento de "tabela_cotas" é uma baliza: baliza[0] é a posição x e baliza[1] é a lista de alturas z.
# Aqui juntamos todas as alturas reais e criamos um grid mais fino, que será usado para
# avaliar as splines das balizas.
# ============================================================

# Lista única com todas as alturas z medidas em todas as balizas
zs_reais = [z for baliza in tabela_cotas for z in baliza[1]]

# Subdivide cada intervalo entre alturas reais em 20 partes, gerando um grid vertical refinado
zs_desejados = gerar_grid_refinado(zs_reais, subdivisoes_por_intervalo=20)

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
# ============================================================

xs_balizas = [baliza[0] for baliza in tabela_cotas]  # posição x de cada baliza original
x_min_global = min(xs_balizas)  # extremidade inicial (menor x)
x_max_global = max(xs_balizas)  # extremidade final (maior x)

# 20 posições x igualmente espaçadas entre as duas extremidades (convertido de array numpy para lista)
xs_desejados = np.linspace(x_min_global, x_max_global, 20).tolist()

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
# Vetores normais e visualização.
# Junta todos os painéis em uma única lista, calcula o vetor normal de cada um
# e gera a visualização 3D interativa em HTML.
# ============================================================

# Lista completa de painéis: boreste + bombordo + fundo + pontas
todos_paineis = paineis + paineis_bombordo + paineis_fundo + paineis_pontas

# Vetor normal de cada painel (necessário para os cálculos hidrostáticos)
paineis_vetor = vetor_normal_painel(todos_paineis)

# Gera o arquivo HTML com o casco em 3D; os painéis de fechamento (fundo e pontas) são agrupados
gerar_visualizacao_3d(
    malha_final,                                          # malha de pontos do casco
    paineis_boreste=paineis,                              # painéis de boreste
    paineis_bombordo=paineis_bombordo,                    # painéis de bombordo
    paineis_fechamento=paineis_fundo + paineis_pontas,    # painéis de fundo, proa e popa
    caminho_saida="visualizacao_casco.html",              # arquivo de saída
)