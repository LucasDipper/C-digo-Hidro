import csv

def _abrir_leitor(caminho_arquivo, arquivo_aberto):
    # utf-8-sig remove o BOM se ele existir (e não faz nada se não existir).
    # csv.Sniffer detecta sozinho se o separador é vírgula ou ponto-e-vírgula.
    amostra = arquivo_aberto.read(2048)
    arquivo_aberto.seek(0)
    try:
        dialeto = csv.Sniffer().sniff(amostra, delimiters=',;')
    except csv.Error:
        dialeto = csv.excel   # padrão: vírgula
    return csv.DictReader(arquivo_aberto, dialect=dialeto)

def ler_balizas(caminho_arquivo):
    balizas_x = {}
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])
            x = float(linha['x_m'])
            balizas_x[indice] = x
    return balizas_x

def ler_linhas_dagua(caminho_arquivo):
    linhas_z = {}
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])
            z = float(linha['z'])
            linhas_z[indice] = z
    return linhas_z

def ler_tabela_cotas(caminho_balizas, caminho_linhas_dagua, caminho_cotas):
    balizas_x = ler_balizas(caminho_balizas)
    linhas_z = ler_linhas_dagua(caminho_linhas_dagua)
    pontos_por_baliza = {}

    with open(caminho_cotas, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_cotas, arquivo)
        for linha in leitor:
            indice_baliza = int(linha['baliza'])
            indice_linha_dagua = int(linha['linha_dagua'])
            y = float(linha['cota_y'])
            z = linhas_z[indice_linha_dagua]
            if indice_baliza not in pontos_por_baliza:
                pontos_por_baliza[indice_baliza] = []
            pontos_por_baliza[indice_baliza].append((z, y))

    tabela_cotas = []
    for indice_baliza in sorted(pontos_por_baliza.keys()):
        x = balizas_x[indice_baliza]
        pontos = sorted(pontos_por_baliza[indice_baliza])
        z_conhecidos = [ponto[0] for ponto in pontos]
        y_conhecidos = [ponto[1] for ponto in pontos]
        tabela_cotas.append((x, z_conhecidos, y_conhecidos))
    return tabela_cotas