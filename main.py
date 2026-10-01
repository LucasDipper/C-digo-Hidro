import numpy as np
from leitura import ler_tabela_cotas
from splines import (construir_linha_baliza, construir_malha_por_contorno,
                        formar_paineis, vetor_normal_painel, 
                        gerar_grid_refinado, espelhar_paineis, fechar_fundo,
                        fechar_pontas, contorno_perfil)
from visualizar import gerar_visualizacao_3d

# --- leitura dos dados reais ---
#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-060-4210\Balizas-4210.csv",
#    r"C:\Users\Vieira\Documents\Casco-060-4210\Linhas_dagua-4210.csv",
#    r"C:\Users\Vieira\Documents\Casco-060-4210\Cotas-4210.csv",
#)

#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-065-4211\Balizas-4211.csv",
#    r"C:\Users\Vieira\Documents\Casco-065-4211\Linhas_dagua-4211.csv",
#    r"C:\Users\Vieira\Documents\Casco-065-4211\Cotas-4211.csv",
#)

#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-070-4212\Balizas-4212.csv",
#    r"C:\Users\Vieira\Documents\Casco-070-4212\Linhas_dagua-4212.csv",
#    r"C:\Users\Vieira\Documents\Casco-070-4212\Cotas-4212.csv",
#)

tabela_cotas = ler_tabela_cotas(
    r"C:\Users\Vieira\Documents\Casco-075-4213\Balizas-4213.csv",
    r"C:\Users\Vieira\Documents\Casco-075-4213\Linhas_dagua-4213.csv",
    r"C:\Users\Vieira\Documents\Casco-075-4213\Cotas-4213.csv",
)

#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\balizas.csv",
#    r"C:\Users\Vieira\Documents\linhas_dagua.csv",
#    r"C:\Users\Vieira\Documents\cotas.csv"
#)

zs_reais = [z for baliza in tabela_cotas for z in baliza[1]]
zs_desejados = gerar_grid_refinado(zs_reais, subdivisoes_por_intervalo=20)

malha_balizas = [construir_linha_baliza(baliza, zs_desejados) for baliza in tabela_cotas]

# Número de pontos ao longo do comprimento em CADA linha d'água. Cada linha d'água
# recebe esses pontos apenas entre a sua própria ponta de popa e ponta de proa,
# então a malha já nasce no contorno real do casco (sem pontos "de enchimento").
N_COLUNAS = 97

# Contorno além das estações FP (x = 0) e AP (x = Lbp), lido do plano de linhas (Fig. 12a/12b):
# lista de (z, x_ponta) em ordem crescente de z, com z em metros. x < 0 fica a vante do FP e
# x > Lbp a ré do AP. Deixe [] para fechar com parede plana na estação (comportamento anterior).
# Valores abaixo são só do 4210W. Proa = leitura aproximada da Fig. 12b (inclinação da roda de proa acima
# do calado; z = 6.502 m é a WL 1.00, onde a roda passa exatamente no FP).
EXT_PROA = [(6.502, 0.0), (8.127, -0.49), (9.753, -0.85)]
EXT_POPA = []     
malha_final = construir_malha_por_contorno(malha_balizas, N_COLUNAS, ext_ini=EXT_PROA, ext_fim=EXT_POPA)

paineis = formar_paineis(malha_final)
paineis_bombordo = espelhar_paineis(paineis)
paineis_fundo = fechar_fundo(malha_final)
paineis_pontas, idx_proa, idx_popa = fechar_pontas(malha_final)

# Pontas abertas: linhas d'água cujo casco continua além da primeira/última estação da tabela
# (y > 0 na coluna extrema). Coluna 0 = proa (FP, x = 0); última coluna = popa (AP).
# Sem o contorno de EXT_PROA / EXT_POPA abaixo, essas linhas terminam numa parede plana na estação.
n_abertas_proa = sum(1 for p in malha_final[0] if p is not None and p[1] > 1e-3)
n_abertas_popa = sum(1 for p in malha_final[-1] if p is not None and p[1] > 1e-3)
if n_abertas_proa:
    print(f"Aviso: {n_abertas_proa} linha(s) d'água da proa terminam cortadas na estação x={malha_final[0][-1][0]:.2f} m (dado da tabela acaba ali).")
if n_abertas_popa:
    print(f"Aviso: {n_abertas_popa} linha(s) d'água da popa terminam cortadas na estação x={malha_final[-1][-1][0]:.2f} m (dado da tabela acaba ali).")

todos_paineis = paineis + paineis_bombordo + paineis_fundo + paineis_pontas
paineis_vetor = vetor_normal_painel(todos_paineis)

# Silhueta lateral (plano de simetria) tirada direto da malha: quilha, proa, convés e popa.
perfil_diametral = {"contorno": contorno_perfil(malha_final)}

gerar_visualizacao_3d(
    malha_final,
    paineis_boreste=paineis,
    paineis_bombordo=paineis_bombordo,
    paineis_fechamento=paineis_fundo + paineis_pontas,
    perfil_diametral=perfil_diametral, 
    caminho_saida="visualizacao_casco.html",
)

#https://pt.scribd.com/document/487650326/Lines-and-offset-of-Series-60
#aaaaaaaaaaa
