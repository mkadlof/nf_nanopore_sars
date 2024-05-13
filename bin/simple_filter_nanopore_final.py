'''
prosty parser puszczamy ivar-a 2 razy - pierwszy raz w opcji strict tak by zapisywal tylko
te ready ktore mapuja sie na priner
drugi raz z opcja -e czyli zapisuje wszystko robi to samo co opcja strict, ale dodatkowo ready
nie mapujace sie na primery tez sa zapisywane ale z poprawionym quality check quality check

nastepnie z opcji strict pobieramy nazwy wlasciwych odczytow i dodajemy do nich ich pary z puszczania z opcji -e

Czego nie ma tutaj to readow ktore choc sa wewnatrz konkretnego amplikonu to ani jeden koniec nie mapuje sie na
zaden z dwoch primerow

Uspojnianie z illumina ? To wprawdzie dziala. ale
1. uzywac dluzszych o 1 SARS1 primers + samtools -a na odczytach interamplikon
2. dla odczyt z fuzji ktorych uzywamy jesli amplikon wypadl to  customowy bed do filtorwania
3.
'''

import pysam
import sys


def read_amplicon_scheme(bed, bed_offset = 1):
    """
    Prosta funkcja to wczytywania bed-a ze schematem primerow i tworzeniem na jego podstawie slownika. Offset determinuje
    o ile bp rozszerzamy amplikon w strone 5' i 3'. Z default offset to 0 i raczej tego nie bede zmienial
    :param bed: sciezka do pliku bed z primerami. Opisane w dokumentacji
    :param bed_offset: int o ile rozszerzamy amplikon. Default 1 zakladamym ze uzytkownik poda poprawnego bed-a a
    my sami rozszerzamy go o 1
    :return: 5  slownikow. Slowniki slownik_amplikonow_with_alt_outer i slownik_amplikonow_with_alt_inner maja jako
    klucz numer amplikonu (1,2,3...). Ktore same sa slownikami z tylko dwoma kluczami (LEFT i RIGHT). Wartosciami tych
    podslownikow jest lista. Najczesciej jedno elementowa zawierajaca granice primeru (zewnetrzna lub wewnetrzna w
    zaleznosci od slownika). Jest lista poniewaz trzymamy tam informacje o ewentualnych primerach alt. tak wiec
    {1}:{LEFT}:[1,3] oznacza ze w slowniku jest amplikon 1 i jego primer left zaczyna sie (lub konczy) na pozycji 1,
    ale jest tez primer alternatywny z poczatkiem (koncem na pozycji 3). Uwaga w slownilu zachowujemy sytuacje z bed-a
    w ktorym zakres jest polotwrty <)
    3 ostatnie slowniki to tylko rzeczy do statystyk gdzie jako klucz mamy numer amplikonu i wartosc zawsze 0 na tym etapie

    Stara pomoc
    # ten slownik jako klucz 3ma nazwe amplikonu jak wyzej np 72
    # ale jako wartosci w odroznieniu do slonika wyzej trzyma kolejny slownik
    # z dwoma kluczami "L" i "R"
    # i dopiero te podlisty 3maja nie wartosc a liste list z grnicami amplikonow w postaci list
    # np dla amplikonu 72 ktory ma dwa primery LEFT jeden w granicach 1,10 i drug i5,20
    # oraz jeden primer prawy w granicach 100-140 mamy taki uklad
    # {'72':{'LEFT': [1, 15], 'RIGHT': [140]}}
    # jak widac dajemy zewnetrzne granice primerow bo
    # ivar w przypadku braku mapowania na oba primery odrzuci ten read gdy nie dajemy flagi e
    # tak wiec chcemy tu uwzglednic ready mapujace się na jeden primer, ktore nie 'dociagnely' do konca


    """
    slownik_amplikonow_with_alt_outer = {} # slowanik z zewnetrznymi granicami amplikonu (czyl 5' primer left i 3' right)
    slownik_amplikonow_with_alt_inner = {} # to samo ale granice wewnetrzne (3' left i 5' right)
    slownik_amplikonow_uzycie = {} # slownik z iloscia odczytow mapujacych sie wewnatrz danego amplikonu
    slownik_amplikonow_uzycie_left = {} # slownik z iloscia odczytow mapujacych sie na primer left danego amplikonu
    slownik_amplikonow_uzycie_right = {} # slownik z iloscia odczytow mapujacych sie na primer right danego amplikonu

    with open(bed) as f:
        for line in f:
            line = line.split()
            # pole 3 moze miec 3 albo 4 elementy bo moga wystepowac amplikony alt
            try:
                _, numer, kierunek, alt = line[3].split('_')
            except:
                _, numer, kierunek = line[3].split('_')

            numer = int(numer)

            if numer not in slownik_amplikonow_with_alt_outer.keys():

                slownik_amplikonow_with_alt_outer[numer] = {} # Ten slownik trzyma zewnetrzne granice amplikonow
                slownik_amplikonow_with_alt_inner[numer] = {} # Ten slownik trzyma wewnetrzne granice amplikonow

                slownik_amplikonow_with_alt_outer[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'] = []

                slownik_amplikonow_with_alt_inner[numer]['LEFT'] = []
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'] = []

                slownik_amplikonow_uzycie[numer] = 0
                slownik_amplikonow_uzycie_left[numer] = 0
                slownik_amplikonow_uzycie_right[numer] = 0

            if 'LEFT' == kierunek :
                slownik_amplikonow_with_alt_outer[numer]['LEFT'].append(int(line[1]) - bed_offset)
                slownik_amplikonow_with_alt_inner[numer]['LEFT'].append(int(line[2]))
            elif 'RIGHT' == kierunek:
                slownik_amplikonow_with_alt_outer[numer]['RIGHT'].append(int(line[2]) + bed_offset)
                slownik_amplikonow_with_alt_inner[numer]['RIGHT'].append(int(line[1]))
            else:
                print('Nie rozpoznany kierunek')

    return slownik_amplikonow_with_alt_outer, slownik_amplikonow_with_alt_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right


def write_reads_strict_inner(initial_bam, final_bam, reject_bam,  statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_with_alt_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, cap = 1000, mapq = 30, length = 300):
    """
    Ta funkcja wyciaga odczyty ktore zaczynaja sie i koncza w w primer. Dodatkowo skoro ta funkcja jest zawsze wykonywana
    na tym etapie wyrzucamy read o jakosci mapowania ponizej 30 i dlugosci ponizej length. To ma na celu przyspieszenie
    obliczen bo takie ready i tak nigdzie nie trafia
    ----------PRIMERL--------------------------PRIMERR----------
    ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-------------
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param cap:
    :return:
    """

    all_reads = pysam.AlignmentFile(initial_bam, "rb")
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_with_alt_outer:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                    begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    if begin_amplikon_L  <= reference_start < end_amplikon_L \
                            and end_amplikon_R  > reference_end >= begin_amplikon_R\
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap and read.mapq >= mapq and (reference_end - reference_start) >= length and read.query_length < (reference_end - reference_start) * 1.1 :
                        slownik_amplikonow_uzycie_left[klucz] += 1
                        slownik_amplikonow_uzycie_right[klucz] += 1
                        slownik_amplikonow_uzycie[klucz] += 1
                        pass_reads.write(read)
                        done = True
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz } both primers \n")
                    elif begin_amplikon_L  <= reference_start < end_amplikon_L \
                            and end_amplikon_R  > reference_end >= begin_amplikon_R\
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] > cap and read.mapq >= mapq and (reference_end - reference_start) >= length:
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz } both primers cap \n")
                        done = True
                    elif read.mapq < mapq:
                        done = True
                        statystyki.write(
                            f"{read.qname}\t{reference_start}\t{reference_end}\tpoor quality\n")
                    elif (reference_end - reference_start) < length:
                        done = True
                        statystyki.write(
                            f"{read.qname}\t{reference_start}\t{reference_end}\ttoo short\n")

        if not done:
        # probowalem read nie trafil do wynikowego skryptu
        # dajemy mu 2-ga szanse w kolejnym pass
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right

def write_reads_overshot(initial_bam, final_bam, reject_bam, statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_with_alt_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, overshoot,  cap):
    """
    Ta funkcja zapisuje ready ktore przestrzelowne sa o 10 w stosunku do pozycji amplikonu (jednego badz obu)
    drugiego primera
    ----------PRIMERL--------------------------PRIMERR--------------
    --------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx----------- (przestrzelone oba primery na zewnatrz)
    --------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx----------------- (przestrzelony primer L na zewnatrz)
    -------------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx--------- (przestrzelony primer R na zewnatrz)
    -------------------xxxxxxxxxxxxxxxxxxxxx------------------------ (przestrzelone primery L i R wewnatrz)
    Granica 10 wynika jest empiryczna ale mozna sie podpudowac ze jest to polowa odlegosci do pozycji ktora jest miedzy
    primerami sasiednich amplikonow (primery roznych amplikonow sa oddalone o ok 40 bp)
    :param initial_bam:
    :param final_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param cap:
    :return:
    """
    all_reads = pysam.AlignmentFile(initial_bam, "rb")
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index

    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_with_alt_outer:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                    begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    if ( begin_amplikon_L - overshoot ) <= reference_start < (end_amplikon_L + overshoot)  \
                            and (end_amplikon_R + overshoot)  > reference_end >= (begin_amplikon_R - overshoot) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        # Aby 'zbalansowac' odczyty prz starcie amplikonu (przez sort-a sa zawsze wczesniej) jesli
                        # primer lewy musi byc uzywany mniej niz cap/2
                        slownik_amplikonow_uzycie_left[klucz] += 1
                        slownik_amplikonow_uzycie_right[klucz] += 1
                        slownik_amplikonow_uzycie[klucz] += 1

                        pass_reads.write(read)
                        done = True
                        statystyki.write(
                            f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} overshoot\n")
                    elif ( begin_amplikon_L - overshoot ) <= reference_start < (end_amplikon_L + overshoot)  \
                            and (end_amplikon_R + overshoot)  > reference_end >= (begin_amplikon_R - overshoot) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] > cap:
                        done = True
                        statystyki.write(
                            f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} overshoot  cap\n")

        if not done:
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()
    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right

def write_reads_fusion_strict(initial_bam, final_bam, statystyki, primer_left_outer, primer_left_inner,
                              primer_middle_left_outer,  primer_middle_left_inner, primer_middle_right_outer,
                              primer_middle_right_inner, primer_right_outer, primer_right_inner, used, case,
                              slownik_amplikonow_uzycie, klucz,  cap, overshoot):
    """
    To jest modyfikacja funkcji write_reads_fusion w ktorej zakladamy ze smieci musza zaczynac sie w primerach
    tzn albo jest fuzja -2 i 0 albo 0 + 2 albo wrecz -2 do + 2 i kazda z tych mozliwosci rozwazamy, dla amplikonow pierwsz
    ego i ostatniego oczywiscie czesc z tych mozliwosci jest zablokowana. Case 1 to przypadek fdzie nie ma primera -2,
    Case 2 gdzie nie ma primera +2 , case 3 gdzie jest primer -2 i +2 . Zeby bylo juz latwiej gdy mamy alty w primerach
    to do tej funkcji dajemy po prosy maksymalne zasiegi danego primera czyli zwykly primer + alt i zmienne primer_left_outer itd
    sa po prostu int-ami
    :param initial_bam:
    :param final_bam:
    :param statystyki:
    :param primer_left_outer:
    :param primer_left_inner:
    :param primer_middle_left_outer:
    :param primer_middle_left_inner:
    :param primer_middle_right_outer:
    :param primer_middle_right_inner:
    :param primer_right_outer:
    :param primer_right_inner:
    :param used:
    :param cap:
    :return:
    """
    my_bam = pysam.AlignmentFile(initial_bam, "rb")
    my_bam_out = pysam.AlignmentFile(final_bam, "wb", template=my_bam)

    amplikon_ilosc = slownik_amplikonow_uzycie[klucz]

    # dla skrajnych amplikonow tych kluczy nie bedzie wiec zrobimy je dummy
    try:
        amplikon_ilosc_left_left = slownik_amplikonow_uzycie[klucz - 2]
    except KeyError:
        amplikon_ilosc_left_left = 0
    try:
        amplikon_ilosc_left = slownik_amplikonow_uzycie[klucz - 1]
    except KeyError:
        amplikon_ilosc_left = 0

    try:
        amplikon_ilosc_right = slownik_amplikonow_uzycie[klucz + 1]
    except KeyError:
        amplikon_ilosc_right = 0

    try:
        amplikon_ilosc_right_right = slownik_amplikonow_uzycie[klucz + 2]
    except KeyError:
        print(f'W slowniku nie ma klucza {klucz + 2}')
        amplikon_ilosc_right_right = 0

    # skoto obejmuje cale amplikony to podbijam tez ich uzycie
    with open(statystyki, 'w') as f:
        # lecimy po parach odczytów
        for read in my_bam:
            if read.qname in used.keys():
                continue
            reference_start = read.reference_start
            reference_end = read.reference_end
            # w odroznieniu od funkcji wyzej tu nie dajemy slownika amplikonow
            # a wybrane listy z zakresami


            if (case == 1 or case == 3) and (primer_middle_left_outer - overshoot ) <= reference_start < (primer_middle_left_inner + overshoot) \
                    and (primer_right_inner - overshoot ) <= reference_end < ( primer_right_outer + overshoot ) and amplikon_ilosc <= cap:
                    # read jest od lewego primera swojego do prawego priemera +2
                my_bam_out.write(read)
                used[read.qname] = ''
                amplikon_ilosc += 1
                amplikon_ilosc_left += 1
                amplikon_ilosc_left_left +=1
                f.write(
                    f"{read.qname}\t{read.reference_start}\t{read.reference_end}\ttwo_amplicons_{klucz}\n")
                continue
            elif (case == 2 or case == 3) and (primer_left_outer - overshoot ) <= reference_start < (primer_left_inner + overshoot ) and \
                    (primer_middle_right_inner - overshoot ) <= reference_end < (primer_middle_right_outer + overshoot ) and amplikon_ilosc <= cap:
                # read jest od lewego primera -2 do swojego primer prawgo
                my_bam_out.write(read)
                used[read.qname] = ''
                amplikon_ilosc += 1
                amplikon_ilosc_right += 1
                amplikon_ilosc_right_right += 1
                f.write(
                    f"{read.qname}\t{read.reference_start}\t{read.reference_end}\ttwo_amplicons_{klucz}\n")
                continue
            elif case == 3 and (primer_left_outer - overshoot) <= reference_start < (primer_left_inner + overshoot ) and \
                    (primer_right_inner - overshoot )<= reference_end < (primer_right_outer + overshoot ) and amplikon_ilosc <= cap:
                # tylko w case 3 pozwalamy na sprawdzanie od -2 do +2
                my_bam_out.write(read)
                used[read.qname] = ''
                amplikon_ilosc += 1
                amplikon_ilosc_left += 1
                amplikon_ilosc_left_left += 1
                amplikon_ilosc_right += 1
                amplikon_ilosc_right_right += 1
                f.write(
                    f"{read.qname}\t{read.reference_start}\t{read.reference_end}\ttwo_amplicons_{klucz}\n")
                continue


    my_bam_out.close()
    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sorted.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sorted.bam")

    # updatujemy slownik uzycia o obecne w nim oryginalne klucze
    if ( klucz - 2) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz - 2 ] = amplikon_ilosc_left_left
    if (klucz - 1) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz - 1] = amplikon_ilosc_left
    if ( klucz) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz] = amplikon_ilosc
    if (klucz + 1) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz + 1] = amplikon_ilosc_right
    if (klucz + 2) in slownik_amplikonow_uzycie.keys():
        slownik_amplikonow_uzycie[klucz + 2] = amplikon_ilosc_right_right

    return slownik_amplikonow_uzycie, used

def write_reads_partstrict_inner(initial_bam, final_bam, reject_bam, statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_with_alt_inner, slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left,
                      slownik_amplikonow_uzycie_right, cap = 1000, overshoot = 0, usage = 0.6):
    """
    Ta funkcja wyciaga odczyty ktorych mapowanie zaczyna w jednym z primerow. Mniej cenne niz write_reads_strict_inner
    ale tez dosc dobrze wiadomo skad pochodza odczyty
    Przywrocilem ta funkcje bo w sampu 9 wyraznie korzystaja z takich readow, gorzej ze nie sa zbyt konsystentni
    DAjemy empitycznie ze read obejmuje jednak co najmniej polowe amplikonu
    ----------PRIMERL--------------------------PRIMERR----------
    ------------xxxxxxxxxxxxxxxxxxxxxxxxxx----------------------
    --------------------xxxxxxxxxxxxxxxxxxxxxxxxxx--------------
    ALA takie odczyty NIE MOGA
    ----------PRIMERL------------------PRIMERLinnyklucz--------------PRIMERR
    ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------------------
    ALE moga byc tak
    --------------xxxxxxxxxxxxxxxxx-----------------------------------------
    po warunkiem ze to ponad 0.5 uzycia maplikonu. Ta funkcja jest strasznie pod EQA

    zrobic czego takiego
    Ba tym etapie filtruje tez ready ktore sa po prosry za krotkiw
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param slownik_amplikonow_uzycie_left:
    :param slownik_amplikonow_uzycie_right:
    :param cap:
    :param length:
    :return:
    """
    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_uzycie:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_L = slownik_amplikonow_with_alt_inner[klucz]['LEFT'][indeks_start]

                    begin_amplikon_R = slownik_amplikonow_with_alt_inner[klucz]['RIGHT'][indeks_end]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    # jesli read zaczyna sie w lewym primer to 'musi' przejsc lewy primer kolejnego amplikonu + offset

                    try:
                        end_amplikon_L_plus = slownik_amplikonow_with_alt_inner[klucz + 1]['LEFT'][0]
                        begin_amplikon_L_plus = slownik_amplikonow_with_alt_outer[klucz + 1]['LEFT'][0]
                    except:
                        end_amplikon_L_plus = 0
                        begin_amplikon_L_plus = 0

                    # jesli read zaczyna sie w primer prawy to musi przejsc pawy primer amplikonu -1 - offset
                    try:
                        begin_amplikon_R_minus =  slownik_amplikonow_with_alt_inner[klucz - 1 ]['RIGHT'][0]
                        end_amplikon_R_minus = slownik_amplikonow_with_alt_outer[klucz - 1]['RIGHT'][0]
                    except:
                        begin_amplikon_R_minus = 10000000
                        end_amplikon_R_minus = 10000000


                    if (begin_amplikon_L - overshoot) <= reference_start < (end_amplikon_L + overshoot )\
                            and begin_amplikon_R >= reference_end \
                            and ( reference_end > end_amplikon_L_plus or  reference_end < begin_amplikon_L_plus) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        # Aby 'zbalansowac' odczyty prz starcie amplikonu (przez sort-a sa zawsze wczesniej) jesli
                        # primer lewy musi byc uzywany mniej niz cap/2
                        amplicon_zakres = set(range(begin_amplikon_L, end_amplikon_R))
                        read_zakres = set(range(reference_start,reference_end))
                        common_zakres = len(read_zakres.intersection(amplicon_zakres)) / float(len(amplicon_zakres))
                        if common_zakres >= usage:
                            slownik_amplikonow_uzycie_left[klucz] += 1
                            slownik_amplikonow_uzycie[klucz] += 1
                            pass_reads.write(read)
                            done = True
                            statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} only left primer\n")
                        else:
                            done = True
                    elif end_amplikon_L < reference_start \
                            and (end_amplikon_R  + overshoot ) > reference_end >= (begin_amplikon_R - overshoot )\
                            and ( reference_start < begin_amplikon_R_minus or  reference_start > end_amplikon_R_minus) \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        amplicon_zakres = set(range(begin_amplikon_L, end_amplikon_R))
                        read_zakres = set(range(reference_start, reference_end))
                        common_zakres = len(read_zakres.intersection(amplicon_zakres)) / float(len(amplicon_zakres))
                        if common_zakres >= usage:
                            slownik_amplikonow_uzycie_right[klucz] += 1
                            slownik_amplikonow_uzycie[klucz] += 1
                            pass_reads.write(read)
                            done = True
                            statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz} only right primer\n")
                        else:
                            done = True
        if not done:
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sorted.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sorted.bam")

    return slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, slownik_amplikonow_uzycie_right


def write_reads_midnight(initial_bam, final_bam, reject_bam, statystyki, slownik_amplikonow_with_alt_outer,
                      slownik_amplikonow_uzycie, cap = 1000):
    """
    Funkcja ktora bierze smiecoi z midnight, tzn odczyty ktores sa wewnatrz amplikonow, ale nie obejmuja primerow,
    to sa odczyty po filtrowaniu na dlugosc
    :param initial_bam:
    :param final_bam:
    :param reject_bam:
    :param statystyki:
    :param slownik_amplikonow_with_alt_outer:
    :param slownik_amplikonow_with_alt_inner:
    :param slownik_amplikonow_uzycie:
    :param cap:
    :return:
    """

    all_reads = pysam.AlignmentFile(initial_bam, "rb", require_index=False)
    pass_reads = pysam.AlignmentFile(final_bam, "wb", template=all_reads)
    reject_reads = pysam.AlignmentFile(reject_bam, "wb", template=all_reads)
    # uwaga pysam bierze koordynaty z bam-a z 1-index na 0-index


    for read in all_reads:
        done = False
        reference_start = read.reference_start
        reference_end = read.reference_end

        for klucz in slownik_amplikonow_with_alt_outer:
            if done:
                continue
            for indeks_start in range(len(slownik_amplikonow_with_alt_outer[klucz]['LEFT'])):
                for indeks_end in range(len(slownik_amplikonow_with_alt_outer[klucz]['RIGHT'])):
                    begin_amplikon_L = slownik_amplikonow_with_alt_outer[klucz]['LEFT'][indeks_start]
                    end_amplikon_R = slownik_amplikonow_with_alt_outer[klucz]['RIGHT'][indeks_end]

                    if begin_amplikon_L  <= reference_start  \
                            and end_amplikon_R  > reference_end \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] <= cap:
                        slownik_amplikonow_uzycie[klucz] += 1
                        pass_reads.write(read)
                        done = True
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz } midnight case \n")
                    elif begin_amplikon_L  <= reference_start  \
                            and end_amplikon_R  > reference_end \
                            and not done \
                            and slownik_amplikonow_uzycie[klucz] > cap:
                        done = True
                        statystyki.write(f"{read.qname}\t{reference_start}\t{reference_end}\tinside_amplicon {klucz } midnight case cap \n")


        if not done:
        # probowalem read nie trafil do wynikowego skryptu
        # dajemy mu 2-ga szanse w kolejnym pass
            reject_reads.write(read)

    pass_reads.close()
    reject_reads.close()

    final_bam_sort = final_bam.split('.')[0]
    pysam.sort("-o", f"{final_bam_sort}_sorted.bam", final_bam)
    pysam.index(f"{final_bam_sort}_sorted.bam")
    return slownik_amplikonow_uzycie



if __name__ == '__main__':
    all_read = sys.argv[1]  # posortowany output z minimap-a
    bed = sys.argv[2]  # plik bed z primerami
    bed_offset = int(sys.argv[3]) # rozszerz bed-a o tyle
    cap = int(sys.argv[4])  # cap na ilosc odczytow mapujacych sie na konkretny amplikon
    length = int(sys.argv[5])  # cap na wielkosc w ten sposob zapobiegam sztucznemu podbijaniu
    mapq = int(sys.argv[6])
    extra_bed_offset = int(sys.argv[7])
    cap_smieci=int(sys.argv[8]) # cap na smieci czyli ready nie obejmujace dwoch primerow
    statystyki = open('Statystyki.txt', 'w')

    #1 Slownik z ampikonami i uzyciami amplkionow
    slownik_amplikonow_with_alt_outer, slownik_amplikonow_with_alt_inner, \
        slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = read_amplicon_scheme(bed = bed, bed_offset = bed_offset)

    midnight = False
    # podmiana cap smiecie na cap w przypadku dlugich amplikonow z midnight, ktore i tak sa ciete na krotsze odcinki przez tagmentaze
    if (slownik_amplikonow_with_alt_outer[1]['RIGHT'][0] - slownik_amplikonow_with_alt_outer[1]['LEFT'][0]) > 600:
        cap_smieci = cap
        midnight = True

    #2 Pierwszy pass odczyty strict
    slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = write_reads_strict_inner(initial_bam =  all_read,
                                                                   final_bam = 'reads_inner_strict.bam',
                                                                   reject_bam = 'reject_first_pass.bam',
                                                                   statystyki=statystyki,
                                                                   slownik_amplikonow_with_alt_outer=slownik_amplikonow_with_alt_outer,
                                                                   slownik_amplikonow_with_alt_inner=slownik_amplikonow_with_alt_inner,
                                                                   slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                                   slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
                                                                   slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
                                                                   mapq=mapq,
                                                                   length=length,
                                                                   cap=cap)
    #print('First pass:')
    #print(slownik_amplikonow_uzycie)
    #print('Uzycia z lewej')
    #print(slownik_amplikonow_uzycie_left)
    #print('Uzycia z prawej')
    #print(slownik_amplikonow_uzycie_right)
    #print('------')
    pysam.index('reads_inner_strict.bam')

    pysam.sort('-o', 'reject_first_pass_sort.bam', 'reject_first_pass.bam')
    pysam.index("reject_first_pass_sort.bam")



    #3 Ready ktore przestrzeliy jeden lub oba primery beda oddzielnie trimowane

    slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = write_reads_overshot(initial_bam='reject_first_pass_sort.bam',
                                                               final_bam='reads_overshot.bam',
                                                               reject_bam='reject_second_pass.bam',
                                                               statystyki=statystyki,
                                                               slownik_amplikonow_with_alt_outer=slownik_amplikonow_with_alt_outer,
                                                               slownik_amplikonow_with_alt_inner=slownik_amplikonow_with_alt_inner,
                                                               slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                               slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
                                                               slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
                                                               overshoot=extra_bed_offset,
                                                               cap=cap)
    #print('Second pass:')
    #print(slownik_amplikonow_uzycie)
    #print('Uzycia z lewej')
    #print(slownik_amplikonow_uzycie_left)
    #print('Uzycia z prawej')
    #print(slownik_amplikonow_uzycie_right)
    #print('------')
    pysam.index('reads_overshot.bam')
    # pysam.index('reject_second_pass.bam')
    pysam.sort('-o', 'reject_second_pass_sort.bam', 'reject_second_pass.bam')
    pysam.index("reject_second_pass_sort.bam")



    #4 Ready smieci nie wiadomo z jakiego sa maplikonu ale mapuja sie w duzej mierze na amplikon o niskim uzyciu
    # moze pochodzi z fuzji sasiednich ampikonow gdy wypada primer (z powodu delecji, albo slabej hybrydyzacji)
    # albo bezposrednio z materialu wyjsciowego generalnie smietnik uzywany do podbicia coverage

    used = {}  # aby zapobiec dublowaniu sie odczytu ajko uzupelnienie ronzych amplikonow musimy to zrobic tak
    # zeby pomoc amplikonom 1 i 2 oraz ostatnim i przedostatnim modyfukjemy slownik_amplikonow_with_alt_outer
    # i slownik_amplikonow_with_alt_inner aby zawrieral primery dla -1 i -2 i + 1 i +2 takze dla nich
    # ale te symulowane granice sa identyczne ze skrajnymi primerami

    min_amplikon = min(slownik_amplikonow_with_alt_outer.keys())
    max_amplikon = max(slownik_amplikonow_with_alt_outer.keys())
    slownik_amplikonow_with_alt_outer[min_amplikon - 1] =   slownik_amplikonow_with_alt_outer[min_amplikon]
    slownik_amplikonow_with_alt_inner[min_amplikon - 1] =   slownik_amplikonow_with_alt_inner[min_amplikon]

    slownik_amplikonow_with_alt_outer[min_amplikon - 2] =   slownik_amplikonow_with_alt_outer[min_amplikon]
    slownik_amplikonow_with_alt_inner[min_amplikon - 2] =   slownik_amplikonow_with_alt_inner[min_amplikon]

    slownik_amplikonow_with_alt_outer[max_amplikon + 1] =   slownik_amplikonow_with_alt_outer[max_amplikon]
    slownik_amplikonow_with_alt_inner[max_amplikon + 1] =   slownik_amplikonow_with_alt_inner[max_amplikon]

    slownik_amplikonow_with_alt_outer[max_amplikon + 2] =   slownik_amplikonow_with_alt_outer[max_amplikon]
    slownik_amplikonow_with_alt_inner[max_amplikon + 2] =   slownik_amplikonow_with_alt_inner[max_amplikon]

    # strict fusion
    # Tez w dwoch krokach najpierw ready 'pewne', potem przestrzelone.
    # W tym wszystkim nie chodzi o kwestie maskowania primer'ow,
    # a o fakt jakiej jakosci ready bierzemy do analizy
    for klucz, wartosc in slownik_amplikonow_uzycie.items():
        if wartosc < cap:
            print(f'Analizuje fuzje ampliku {klucz} z sasiednimi amplikonami')
            initial_bam = "reject_second_pass_sort.bam"
            final_bam = f'reads_two_amplicons_{klucz}.bam'
            statystyki_new = f'Statystyki_two_amplicons_{klucz}.txt'
            # slownik uzycia ma normalnie numerowane primery, a slowniki slownik_amplikonow_with_alt_outer
            # slownik_amplikonow_with_alt_inner i inner maja dummy o 2 rozszerzone wiec ponizej zawsze dziala


            slownik_amplikonow_uzycie, used = write_reads_fusion_strict(initial_bam=initial_bam,
                                                                        final_bam=final_bam,
                                                                        statystyki=statystyki_new,
                                                                        primer_left_outer=min(slownik_amplikonow_with_alt_outer[(klucz - 2)]['LEFT']),
                                                                        primer_left_inner=max(slownik_amplikonow_with_alt_inner[(klucz - 2)]['LEFT']),
                                                                        primer_middle_left_outer=min(slownik_amplikonow_with_alt_outer[(klucz)]['LEFT']),
                                                                        primer_middle_left_inner=max(slownik_amplikonow_with_alt_inner[(klucz)]['LEFT']),
                                                                        primer_middle_right_outer=max(slownik_amplikonow_with_alt_outer[(klucz)]['RIGHT']),
                                                                        primer_middle_right_inner=min(slownik_amplikonow_with_alt_inner[(klucz)]['RIGHT']),
                                                                        primer_right_outer=max(slownik_amplikonow_with_alt_outer[(klucz + 2)]['RIGHT']),
                                                                        primer_right_inner=min(slownik_amplikonow_with_alt_inner[(klucz + 2)]['RIGHT']),
                                                                        slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                                        klucz=klucz,
                                                                        used=used,
                                                                        case=3,
                                                                        overshoot = extra_bed_offset,
                                                                        cap=cap)


    #print('Third pass:')
    #print(slownik_amplikonow_uzycie)
    #print('Uzycia z lewej')
    #print(slownik_amplikonow_uzycie_left)
    #print('Uzycia z prawej')
    #print(slownik_amplikonow_uzycie_right)




# Przywracamy ready z czesiowym mapowaniem
# subroutine na midnight

    if midnight:
        usage = 0.2
        pysam.collate('-o', "reject_second_pass_sort_unsort.bam", 'reject_second_pass_sort.bam')
        initial_bam = 'reject_second_pass_sort_unsort.bam'
    else:
        initial_bam = 'reject_second_pass_sort.bam'
        usage = 0.5

    slownik_amplikonow_uzycie, slownik_amplikonow_uzycie_left, \
        slownik_amplikonow_uzycie_right = write_reads_partstrict_inner(initial_bam=initial_bam,
                                                                       final_bam='reads_partial_strict.bam',
                                                                       reject_bam='reject_partial_strict_inner.bam',
                                                                       statystyki=statystyki,
                                                                       slownik_amplikonow_with_alt_outer=slownik_amplikonow_with_alt_outer,
                                                                       slownik_amplikonow_with_alt_inner=slownik_amplikonow_with_alt_inner,
                                                                       slownik_amplikonow_uzycie=slownik_amplikonow_uzycie,
                                                                       slownik_amplikonow_uzycie_left=slownik_amplikonow_uzycie_left,
                                                                       slownik_amplikonow_uzycie_right=slownik_amplikonow_uzycie_right,
                                                                       overshoot=extra_bed_offset,
                                                                       usage = usage,
                                                                       cap=cap_smieci)


    if midnight:
        pysam.collate('-o', "reject_partial_strict_inner_unsort.bam", 'reject_partial_strict_inner.bam')
        #pysam.sort('-o', 'reject_partial_strict_inner_sort.bam', 'reject_partial_strict_inner.bam')
        #pysam.index("reject_partial_strict_inner_sort.bam")

        slownik_amplikonow_uzycie =  write_reads_midnight(initial_bam = 'reject_partial_strict_inner_unsort.bam',
                                                          final_bam = 'reads_smieci.bam',
                                                          reject_bam = 'final_refect.bam',
                                                          statystyki = statystyki,
                                                          slownik_amplikonow_with_alt_outer = slownik_amplikonow_with_alt_outer,
                                                          slownik_amplikonow_uzycie = slownik_amplikonow_uzycie,
                                                          cap=cap_smieci)

    with open('output_primers_summary.txt', 'w') as f:
        f.write('Amplikon ID\tAmplion Usage\tPrimers Usage\n')
        for klucz, wartosc in slownik_amplikonow_uzycie.items():
            amplion_usage = slownik_amplikonow_uzycie_left[klucz] + slownik_amplikonow_uzycie_right[klucz]
            f.write(f'{klucz}\t{wartosc}\t{amplion_usage}\n')
    #print('Fourth pass:')
    #print(slownik_amplikonow_uzycie)
    #print('Uzycia z lewej')
    #print(slownik_amplikonow_uzycie_left)
    #print('Uzycia z prawej')
    #print(slownik_amplikonow_uzycie_right)
    #print('------')
    #pysam.index('reads_partial_strict.bam')
    # pysam.index('reject_second_pass.bam')
    pysam.sort('-o', 'reads_partial_strict_sort.bam', 'reads_partial_strict.bam')
    pysam.index("reads_partial_strict_sort.bam")




