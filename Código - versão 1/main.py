import numpy as np
from leitura import ler_tabela_cotas
from splines import (construir_linha_baliza, construir_linha_dagua,
                        montar_malha_final, formar_paineis, vetor_normal_painel, 
                        gerar_grid_refinado, espelhar_paineis, fechar_fundo,
                        fechar_pontas, indice_coluna_valida)
from visualizar import gerar_visualizacao_3d

# --- leitura dos dados reais ---
tabela_cotas = ler_tabela_cotas(
    r"C:\Users\Vieira\Documents\Casco-060-4210\Balizas-4210.csv",
    r"C:\Users\Vieira\Documents\Casco-060-4210\Linhas_dagua-4210.csv",
    r"C:\Users\Vieira\Documents\Casco-060-4210\Cotas-4210.csv",
)

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

#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\Casco-075-4213\Balizas-4213.csv",
#    r"C:\Users\Vieira\Documents\Casco-075-4213\Linhas_dagua-4213.csv",
#    r"C:\Users\Vieira\Documents\Casco-075-4213\Cotas-4213.csv",
#)

#tabela_cotas = ler_tabela_cotas(
#    r"C:\Users\Vieira\Documents\balizas.csv",
#    r"C:\Users\Vieira\Documents\linhas_dagua.csv",
#    r"C:\Users\Vieira\Documents\cotas.csv"
#)

zs_reais = [z for baliza in tabela_cotas for z in baliza[1]]
zs_desejados = gerar_grid_refinado(zs_reais, subdivisoes_por_intervalo=20)

malha_balizas = [construir_linha_baliza(baliza, zs_desejados) for baliza in tabela_cotas]

xs_balizas = [baliza[0] for baliza in tabela_cotas]
x_min_global = min(xs_balizas)
x_max_global = max(xs_balizas)
xs_desejados = np.linspace(x_min_global, x_max_global, 20).tolist()

malha_linhas_dagua = [construir_linha_dagua(malha_balizas, j, xs_desejados)
                       for j in range(len(zs_desejados))]

malha_final = montar_malha_final(malha_linhas_dagua, len(xs_desejados))

paineis = formar_paineis(malha_final)
paineis_bombordo = espelhar_paineis(paineis)
paineis_fundo = fechar_fundo(malha_final)
paineis_pontas, idx_proa, idx_popa = fechar_pontas(malha_final)

if idx_popa != len(malha_final) - 1:
    print(f"Aviso: dado insuficiente na popa. Fechamento feito na baliza x={xs_desejados[idx_popa]:.2f} m, "
          f"em vez de x={xs_desejados[-1]:.2f} m (faltam {xs_desejados[-1]-xs_desejados[idx_popa]:.2f} m sem representar).")

if idx_proa != 0:
    print(f"Aviso: dado insuficiente na proa. Fechamento feito na baliza x={xs_desejados[idx_proa]:.2f} m.")

todos_paineis = paineis + paineis_bombordo + paineis_fundo + paineis_pontas
paineis_vetor = vetor_normal_painel(todos_paineis)

gerar_visualizacao_3d(
    malha_final,
    paineis_boreste=paineis,
    paineis_bombordo=paineis_bombordo,
    paineis_fechamento=paineis_fundo + paineis_pontas,
    caminho_saida="visualizacao_casco.html",
)
#aaaaaaaaa