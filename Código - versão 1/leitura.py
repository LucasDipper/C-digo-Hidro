import csv   # Biblioteca padrão do Python para ler arquivos CSV

# ============================================================
# LEITURA DE ARQUIVOS CSV
# Este arquivo lê os três CSVs de cada embarcação (balizas, linhas d'água e cotas) e
# monta a "tabela_cotas", que é a entrada de construir_linha_baliza() em splines.py.
# Colunas esperadas em cada CSV:
#   - balizas:      "indice" e "x_m"             (posição x de cada baliza, em metros)
#   - linhas d'água: "indice" e "z"              (altura z de cada linha d'água)
#   - cotas:        "baliza", "linha_dagua" e "cota_y"   (y medido em cada cruzamento baliza x linha d'água)
# ============================================================

# ============================================================
# FUNÇÃO AUXILIAR: cria o leitor de CSV
# Usada pelas três funções de leitura, para tratar do mesmo jeito arquivos com
# separador vírgula ou ponto-e-vírgula.
# ============================================================

def _abrir_leitor(caminho_arquivo, arquivo_aberto):
    # Recebe um arquivo JÁ aberto e devolve um leitor (DictReader) que entrega cada linha como dicionário
    # O encoding utf-8-sig, que remove o BOM (marca invisível no início do arquivo, comum em CSVs
    # salvos pelo Excel) se ele existir, é definido no open() de quem chama esta função.
    # csv.Sniffer detecta sozinho se o separador é vírgula ou ponto-e-vírgula.
    amostra = arquivo_aberto.read(2048)   # lê só o começo do arquivo para o Sniffer analisar
    arquivo_aberto.seek(0)                # volta ao início, para o DictReader ler o arquivo todo
    try:
        dialeto = csv.Sniffer().sniff(amostra, delimiters=',;')   # tenta descobrir o separador (',' ou ';')
    except csv.Error:
        dialeto = csv.excel   # padrão: vírgula (usado quando o Sniffer não consegue decidir)
    return csv.DictReader(arquivo_aberto, dialect=dialeto)   # a 1ª linha do CSV vira o nome das colunas

# ============================================================
# LEITURA DAS BALIZAS E DAS LINHAS D'ÁGUA
# Cada uma devolve um dicionário {índice: valor}, para depois ligar os índices da tabela de cotas
# às posições reais (x das balizas, z das linhas d'água).
# ============================================================

def ler_balizas(caminho_arquivo):
    # Lê o CSV das balizas e devolve {índice da baliza: posição x em metros}.
    balizas_x = {}   # dicionário de resultado
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:   # "with" fecha o arquivo sozinho no fim
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:                  # cada linha é um dicionário com as colunas do CSV
            indice = int(linha['indice'])     # número da baliza (inteiro)
            x = float(linha['x_m'])           # posição longitudinal da baliza, em metros
            balizas_x[indice] = x
    return balizas_x

def ler_linhas_dagua(caminho_arquivo):
    # Lê o CSV das linhas d'água e devolve {índice da linha d'água: altura z}.
    linhas_z = {}    # dicionário de resultado
    with open(caminho_arquivo, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_arquivo, arquivo)
        for linha in leitor:
            indice = int(linha['indice'])     # número da linha d'água (inteiro)
            z = float(linha['z'])             # altura da linha d'água
            linhas_z[indice] = z
    return linhas_z

# ============================================================
# MONTAGEM DA TABELA DE COTAS
# Combina os três arquivos. Para cada baliza, reúne os pares (z, y) medidos e devolve
# uma lista de tuplas (x, lista de z, lista de y), uma por baliza, no formato esperado por
# construir_linha_baliza() em splines.py.
# ============================================================

def ler_tabela_cotas(caminho_balizas, caminho_linhas_dagua, caminho_cotas):
    balizas_x = ler_balizas(caminho_balizas)               # {índice da baliza: x}
    linhas_z = ler_linhas_dagua(caminho_linhas_dagua)      # {índice da linha d'água: z}
    pontos_por_baliza = {}   # {índice da baliza: lista de pares (z, y) medidos nessa baliza}

    # --- LEITURA DAS COTAS ---
    with open(caminho_cotas, newline='', encoding='utf-8-sig') as arquivo:
        leitor = _abrir_leitor(caminho_cotas, arquivo)
        for linha in leitor:
            indice_baliza = int(linha['baliza'])                 # a qual baliza pertence esta cota
            indice_linha_dagua = int(linha['linha_dagua'])       # em qual linha d'água foi medida
            y = float(linha['cota_y'])                           # valor da cota y (meia-boca)
            z = linhas_z[indice_linha_dagua]                     # converte o índice da linha d'água na altura z real
            if indice_baliza not in pontos_por_baliza:
                pontos_por_baliza[indice_baliza] = []            # primeira cota desta baliza: cria a lista
            pontos_por_baliza[indice_baliza].append((z, y))      # guarda o par (z, y)

    # --- MONTAGEM DA TABELA FINAL ---
    tabela_cotas = []   # lista de tuplas (x, z_conhecidos, y_conhecidos), uma por baliza
    for indice_baliza in sorted(pontos_por_baliza.keys()):                # balizas em ordem crescente de índice
        x = balizas_x[indice_baliza]                                      # posição x da baliza
        pontos = sorted(pontos_por_baliza[indice_baliza])                 # ordena os pares por z crescente (o que splines.py espera)
        z_conhecidos = [ponto[0] for ponto in pontos]                     # separa as alturas z
        y_conhecidos = [ponto[1] for ponto in pontos]                     # separa as cotas y correspondentes
        tabela_cotas.append((x, z_conhecidos, y_conhecidos))
    return tabela_cotas