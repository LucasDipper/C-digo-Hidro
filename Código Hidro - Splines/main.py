import numpy as np  # Biblioteca numérica (arrays, linspace etc.)
from leitura import ler_tabela_cotas  # Função que lê os 3 CSVs e monta a tabela de cotas
from splines import (construir_linha_baliza, construir_malha_por_contorno,
                        formar_paineis, vetor_normal_painel, 
                        gerar_grid_refinado, espelhar_paineis, fechar_fundo,
                        fechar_pontas, contorno_perfil)  # Funções de interpolação (splines), malha, painéis e perfil
from visualizar import gerar_visualizacao_3d  # Função que gera o HTML com a visualização 3D do casco
#A parte acima importa bibliotecas do python e funções dos outros código/arquivos dessa pasta

# ============================================================
# Leitura de arquivos, balizas, linhas d'agua e cotas das embarcações selecionadas.
# É chamado a função "ler_tabela_cotas()", os dados são enviados para a variavel "tabela_cotas".
# Só uma embarcação deve estar ativa por vez: as demais ficam comentadas (#).
# Para trocar de casco, comente o bloco ativo e descomente o desejado.
# ============================================================

# --- leitura dos dados reais ---
# Argumentos de cada chamada: (arquivo das balizas, arquivo das linhas d'água, arquivo das cotas).
# O prefixo r"..." indica "raw string", para o Python não interpretar as barras invertidas do caminho do Windows.

#Embarcação série 60. Coeficiente de bloco: 0.60
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
# Refinamento do grid vertical (alturas z) e splines das balizas.
# Cada elemento de "tabela_cotas" é uma baliza: baliza[0] é a posição x e baliza[1] é a lista de alturas z.
# Juntamos todas as alturas reais, criamos um grid mais fino e, para cada baliza,
# ajustamos uma spline y(z) avaliada nesse grid, obtendo uma curva suave por seção.
# ============================================================

# Lista única com todas as alturas z medidas em todas as balizas (list comprehension aninhada)
zs_reais = [z for baliza in tabela_cotas for z in baliza[1]]

# Subdivide cada intervalo entre alturas reais em 20 partes, gerando um grid vertical refinado
zs_desejados = gerar_grid_refinado(zs_reais, subdivisoes_por_intervalo=20)

# Uma linha de baliza (curva suave) para cada baliza da tabela, avaliada em zs_desejados
malha_balizas = [construir_linha_baliza(baliza, zs_desejados) for baliza in tabela_cotas]

# ============================================================
# Malha por contorno (parte nova desta versão).
# Em vez de usar o mesmo grid de posições x para todas as alturas, cada linha d'água recebe
# seus pontos apenas entre a sua própria ponta de popa e ponta de proa. Assim, a malha já
# nasce no contorno real do casco, e o contorno da proa/popa além das estações
# extremas (FP e AP) pode ser informado por listas (EXT_PROA e EXT_POPA).
# ============================================================

# Número de pontos ao longo do comprimento em CADA linha d'água. Cada linha d'água
# recebe esses pontos apenas entre a sua própria ponta de popa e ponta de proa,
# então a malha já nasce no contorno real do casco (sem pontos "de enchimento").
N_COLUNAS = 97

# Contorno além das estações FP (x = 0) e AP (x = Lbp), lido do plano de linhas (Fig. 12a/12b):
# lista de (z, x_ponta) em ordem crescente de z, com z em metros. x < 0 fica a vante do FP e
# x > Lbp a ré do AP. Deixe [] para fechar com parede plana na estação (comportamento anterior).
# Valores abaixo são só do 4210W. Proa = leitura aproximada da Fig. 12b (inclinação da roda de proa acima
# do calado; z = 6.502 m é a WL 1.00, onde a roda passa exatamente no FP).
# Atenção: se trocar de casco acima, revise estes valores (ou deixe []), pois não valem para os outros cascos.
EXT_PROA = [(6.502, 0.0), (8.127, -0.49), (9.753, -0.85)]   # proa: pares (z, x_ponta); x negativo = à frente do FP
EXT_POPA = []     # popa: sem contorno informado, fecha com parede plana na estação AP

# Constrói a malha final: para cada linha d'água, N_COLUNAS pontos entre as pontas de popa e proa,
# usando as balizas ajustadas e os contornos EXT_PROA / EXT_POPA além das estações extremas.
malha_final = construir_malha_por_contorno(malha_balizas, N_COLUNAS, ext_ini=EXT_PROA, ext_fim=EXT_POPA)

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

# ============================================================
# Verificação das pontas abertas.
# Verifica se alguma linha d'água termina "cortada" na estação extrema, ou seja, o casco
# continua além do último dado da tabela, mas a malha acaba ali.
# ============================================================

# Pontas abertas: linhas d'água cujo casco continua além da primeira/última estação da tabela
# (y > 0 na coluna extrema). Coluna 0 = proa (FP, x = 0); última coluna = popa (AP).
# Sem o contorno de EXT_PROA / EXT_POPA acima, essas linhas terminam numa parede plana na estação.
# Conta quantos pontos da coluna extrema existem (não são None) e têm y acima de 1e-3 (tolerância para "zero"):
n_abertas_proa = sum(1 for p in malha_final[0] if p is not None and p[1] > 1e-3)    # primeira coluna (proa)
n_abertas_popa = sum(1 for p in malha_final[-1] if p is not None and p[1] > 1e-3)   # última coluna (popa)
# Se houver linhas abertas, avisa quantas são e em qual x a malha termina (x do ponto mais alto da coluna)
if n_abertas_proa:
    print(f"Aviso: {n_abertas_proa} linha(s) d'água da proa terminam cortadas na estação x={malha_final[0][-1][0]:.2f} m (dado da tabela acaba ali).")
if n_abertas_popa:
    print(f"Aviso: {n_abertas_popa} linha(s) d'água da popa terminam cortadas na estação x={malha_final[-1][-1][0]:.2f} m (dado da tabela acaba ali).")

# ============================================================
# Vetores normais, perfil e visualização.
# Junta todos os painéis em uma única lista, calcula o vetor normal de cada um,
# extrai a silhueta lateral do casco e gera a visualização 3D interativa em HTML.
# ============================================================

# Lista completa de painéis: boreste + bombordo + fundo + pontas
todos_paineis = paineis + paineis_bombordo + paineis_fundo + paineis_pontas

# Vetor normal de cada painel (necessário para os cálculos hidrostáticos)
paineis_vetor = vetor_normal_painel(todos_paineis)

# Silhueta lateral (plano de simetria) tirada direto da malha: quilha, proa, convés e popa.
# É guardada em um dicionário com a chave "contorno", que o visualizador usa para desenhar o perfil.
perfil_diametral = {"contorno": contorno_perfil(malha_final)}

# Gera o arquivo HTML com o casco em 3D; os painéis de fechamento (fundo e pontas) são agrupados
gerar_visualizacao_3d(
    malha_final,                                          # malha de pontos do casco
    paineis_boreste=paineis,                              # painéis de boreste
    paineis_bombordo=paineis_bombordo,                    # painéis de bombordo
    paineis_fechamento=paineis_fundo + paineis_pontas,    # painéis de fundo, proa e popa
    perfil_diametral=perfil_diametral,                    # silhueta lateral do casco (novo nesta versão)
    caminho_saida="visualizacao_casco.html",              # arquivo de saída
)

#https://pt.scribd.com/document/487650326/Lines-and-offset-of-Series-60
#aaaaaaaaaaa
