import csv   # módulo padrão do Python para ler arquivos CSV - não precisa instalar nada

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
    # Lê balizas.csv e devolve um dicionário {indice: x_m},
    # para depois conseguir traduzir "baliza número tal" em "posição x tal".
    balizas_x = {}
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        # csv.DictReader lê cada linha como um dicionário, usando o cabeçalho
        # (a primeira linha do CSV) como as chaves - assim você acessa por
        # nome de coluna (linha['x_m']) em vez de por posição (linha[2])
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])   # tudo que vem do CSV chega como texto - precisa converter
            x = float(linha['x_m'])
            balizas_x[indice] = x
    return balizas_x
    # o "with open(...) as arquivo" garante que o arquivo é fechado sozinho
    # no final, mesmo se der erro no meio - não precisa chamar arquivo.close()

def ler_linhas_dagua(caminho_arquivo):
    # Mesma lógica de ler_balizas, mas para linhas_dagua.csv -> {indice: z}
    linhas_z = {}
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])
            z = float(linha['z'])
            linhas_z[indice] = z
    return linhas_z

def ler_tabela_cotas(caminho_balizas, caminho_linhas_dagua, caminho_cotas):
    # Junta os três arquivos e monta a tabela no formato que o resto do
    # código já entende: uma lista de tuplas (x_baliza, [z conhecidos], [y conhecidos]).
    balizas_x = ler_balizas(caminho_balizas)
    linhas_z = ler_linhas_dagua(caminho_linhas_dagua)

    # agrupa os pontos (z, y) de cada baliza num dicionário temporário:
    # {indice_baliza: [(z1, y1), (z2, y2), ...]}
    pontos_por_baliza = {}

    with open(caminho_cotas, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_cotas, arquivo)
        for linha in leitor:
            indice_baliza = int(linha['baliza'])
            indice_linha_dagua = int(linha['linha_dagua'])
            y = float(linha['cota_y'])

            z = linhas_z[indice_linha_dagua]   # traduz o índice da linha d'água para a altura real em metros

            if indice_baliza not in pontos_por_baliza:
                pontos_por_baliza[indice_baliza] = []   # primeira vez que essa baliza aparece: cria a lista
            pontos_por_baliza[indice_baliza].append((z, y))

    # monta a lista final, uma tupla por baliza, com os pontos ordenados por z crescente
    # (necessário porque cotas.csv nem sempre lista as linhas d'água em ordem, e a spline exige isso)
    tabela_cotas = []
    for indice_baliza in sorted(pontos_por_baliza.keys()):
        x = balizas_x[indice_baliza]
        pontos = sorted(pontos_por_baliza[indice_baliza])   # ordena as tuplas (z,y) pelo z (primeiro item)
        z_conhecidos = [ponto[0] for ponto in pontos]
        y_conhecidos = [ponto[1] for ponto in pontos]
        tabela_cotas.append((x, z_conhecidos, y_conhecidos))

    return tabela_cotas
