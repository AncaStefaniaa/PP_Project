#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

struct px {
    unsigned char b, g, r;
};

struct fi {
    unsigned int x, y, W, H, culoare;
    bool ok;
    double correlation;

};

void incarcareImg2(char *bmpIn, unsigned int *H, unsigned int *W, struct px ***liniarizare, unsigned char **header) {

    FILE *fin = fopen(bmpIn, "rb");
    fin = fopen(bmpIn, "rb");

    if (fin == NULL) {
        printf("nu am gasit imaginea sursa din care citesc");
        return;
    }

    //citesc dimensiunile de la octetul 18 respectiv 22
    fseek(fin, 18, SEEK_SET);
    fread(W, sizeof(unsigned int), 1, fin);
    fread(H, sizeof(unsigned int), 1, fin);
    printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n", (*W), (*H));
    rewind(fin);

    int padding;
    *header = (unsigned char *) malloc(54 * sizeof(unsigned char));
    *liniarizare = (struct px **) malloc((*H) * sizeof(struct px *));

    fread(*header, 1, 54, fin);

    fseek(fin, 54, SEEK_SET);
    int i, j, ct = 0;
    padding = (4 - (*W * 3) % 4) % 4;
    for (i = *H - 1; i >= 0; i--) {
        *(*liniarizare + i) = (struct px *) malloc((*W) * sizeof(struct px));
        for (j = 0; j < *W; j++) {
            fread(&(*liniarizare)[i][j], 1, 3, fin);
            //printf("%d\n",(*liniarizare)[i][j].r);
        }

        fseek(fin, padding * sizeof(char), SEEK_CUR);
    }

    fclose(fin);

}

void grayscale(char *nume_fisier_sursa, unsigned int H, unsigned int W, struct px ***liniarizare) {

    int i, j;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            unsigned char aux;
            aux = 0.299 * ((*liniarizare)[i][j]).r + 0.587 * ((*liniarizare)[i][j]).g +
                  0.114 * ((*liniarizare)[i][j]).b;
            ((*liniarizare)[i][j]).r = ((*liniarizare)[i][j]).g = ((*liniarizare)[i][j]).b = aux;

        }

}


void writeBmpOut2(struct px **liniarizare, char *bmpOut1, char *bmpIn, unsigned char *header, unsigned int latime_img,
                  unsigned int inaltime_img) {

    int i, j;
    FILE *fout = fopen(bmpOut1, "wb");
    rewind(fout);
    //scriu headerul in fisierul binar
    fwrite(header, 1, 54, fout);

    int padding = (4 - (latime_img * 3) % 4) % 4;
    for (i = inaltime_img - 1; i >= 0; i--) {

        for (j = 0; j < latime_img; j++) {
            fwrite(&liniarizare[i][j], 1, 3, fout);
        }

        fwrite((&liniarizare[i][j]) - 1, sizeof(char), padding, fout);
    }

    fclose(fout);
}

void corelatie(struct px **liniarizare, struct px **cifra, double ps, struct fi **D, struct fi *fiCurent,
               unsigned int H, unsigned int W, unsigned int H1, unsigned int W1, unsigned int *numarDetectii,unsigned int *q) {

    unsigned int i, j, x, y, n = 165;
    unsigned char Fi, S;
    double SMed, FiMediu, sumaS, sumaFi;


    for (x = 0; x < H - H1; x++) {
        for (y = 0; y < W - W1; y++) {

            double corr = 0;
            SMed = 0;
            FiMediu = 0;

            for (i = x; i < x + H1; i++)
                for (j = y; j < y + W1; j++) {
                    //printf("%d %d %d %dmesaj\n", i, j, x, y);
                    Fi = liniarizare[i][j].r;
                    S = cifra[i - x][j - y].r;

                    FiMediu += Fi;
                    SMed += S;
                }
            SMed = SMed / n;
            FiMediu = FiMediu / n;

            double devStandardS = 0;
            double devStandardFi = 0;

            for (i = x; i < x + H1; i++)
                for (j = y; j < y + W1; j++) {
                    Fi = liniarizare[i][j].r;
                    S = cifra[i - x][j - y].r;
                    devStandardFi += pow((Fi - FiMediu), 2);
                    devStandardS += pow((S - SMed), 2);
                }
            devStandardFi = sqrt(devStandardFi / (n - 1));
            devStandardS = sqrt(devStandardS / (n - 1));

            for (i = x; i < x + H1; i++)
                for (j = y; j < y + W1; j++) {
                    Fi = liniarizare[i][j].r;
                    S = cifra[i - x][j - y].r;

                    corr += ((Fi - FiMediu) * (S - SMed)) / (devStandardFi * devStandardS);

                }
            corr /= n;


            if (corr >= ps) {
                (*fiCurent).x = x;
                (*fiCurent).y = y;
                (*fiCurent).correlation = corr;
                (*fiCurent).W = W1;
                (*fiCurent).H = H1;
                (*fiCurent).ok = true;
                (*fiCurent).culoare = *q;
                (*numarDetectii)++;
                //printf("%d\n", *numarDetectii);
                //printf("%d\n",(*fiCurent).culoare);
                //printf ("%d\n", *numarDetectii);
                (*D)[*numarDetectii - 1] = (*fiCurent);

            }

        }
    }
/*
    printf("\n---Printam corelatiile detectiilor ---\n");

    for (i = 1; i <= *numarDetectii; i++)
        printf("%.2f; ", (*D)[i].correlation);
*/
    (*q)++;//trecem la urmatoarea culoare

}

void eliminareMaxime(struct fi *D, struct px ***liniarizare, struct px *culoare, unsigned int numarDetectii, unsigned int H1,
                     unsigned int W1,unsigned int *contorSuprapuneri,struct fi **Dsuprapuneri) {

    int i, j, R1Top, R2Top, R1Bootom, R2Bottom;
    for (i = 0; i < numarDetectii - 1; i++)
        for (j = i + 1; j < numarDetectii; j++) {

            double suprapunere = 0;
            int surface1, surface2;
            surface1 = surface2 = H1 * W1;
            double surfacesup = 0;

            //pastrez coordonate detectiei cu corelatia cea mai mare
            R1Top = D[i].x > D[j].x ? D[i].x : D[j].x;
            R2Top = D[i].y > D[j].y ? D[i].y + 11 : D[j].y;
            R1Bootom = D[i].x + 15 < D[j].x + 15 ? D[i].x + 15 : D[j].x + 15;
            R2Bottom = D[i].y + 11 < D[j].y + 11 ? D[i].y + 11 : D[j].y + 11;

            int dif1 = R1Bootom - R1Top;
            int dif2 = R2Bottom - R2Top;
            surfacesup = (double) ((dif1 + 1) * (dif2 + 1)) / 3;
            if ((R1Bootom < R1Top) || (R2Bottom < R2Top))
                surfacesup = 0;

            if (surface1 + surface2 - surfacesup)
                suprapunere = surfacesup / (surface1 + surface2 - surfacesup);
            if (suprapunere > 0.2)
                (*Dsuprapuneri)[*contorSuprapuneri++]=D[j];


        }

}

void initializareCulori(struct px **culoare) {
    int i;
    *culoare = (struct px *) malloc(10 * sizeof(struct px));

    //rosu
    (*culoare)[0].r = 255;
    (*culoare)[0].g = 0;
    (*culoare)[0].b = 0;
    //galben
    (*culoare)[1].r = 255;
    (*culoare)[1].g = 255;
    (*culoare)[1].b = 0;
    //verde
    (*culoare)[2].r = 0;
    (*culoare)[2].g = 255;
    (*culoare)[2].b = 0;
    //cyan
    (*culoare)[3].r = 0;
    (*culoare)[3].g = 255;
    (*culoare)[3].b = 255;
    //magenta
    (*culoare)[4].r = 255;
    (*culoare)[4].g = 0;
    (*culoare)[4].b = 255;
    //albastru
    (*culoare)[5].r = 0;
    (*culoare)[5].g = 0;
    (*culoare)[5].b = 255;
    //argintiu
    (*culoare)[6].r = 192;
    (*culoare)[6].g = 192;
    (*culoare)[6].b = 192;
    //albastru2
    (*culoare)[7].r = 255;
    (*culoare)[7].g = 140;
    (*culoare)[7].b = 0;
    //magenta2
    (*culoare)[8].r = 128;
    (*culoare)[8].g = 0;
    (*culoare)[8].b = 128;
    //albastru3
    (*culoare)[9].r = 128;
    (*culoare)[9].g = 0;
    (*culoare)[9].b = 0;

}

void desenareDreptunghiuri(struct px ***liniarizare, struct fi D, struct px culoare, unsigned int numarDetectii) {

    int i, j;

    int e = D.y;
    int w = D.x;

    printf("Coordonatele laturilor: %d %d %d %d\n", w, w + 15, e, e + 11);

    for (j = D.y; j <= 11 + D.y; j++) {
        (*liniarizare)[w][j].r = culoare.r;
        (*liniarizare)[w][j].g = culoare.g;
        (*liniarizare)[w][j].b = culoare.b;

        (*liniarizare)[w + 15][j].r = culoare.r;
        (*liniarizare)[w + 15][j].g = culoare.g;
        (*liniarizare)[w + 15][j].b = culoare.b;
    }

    for (i = D.x; i <= 15 + D.x; i++) {

        (*liniarizare)[i][e].r = culoare.r;
        (*liniarizare)[i][e].g = culoare.g;
        (*liniarizare)[i][e].b = culoare.b;

        (*liniarizare)[i][e + 11].r = culoare.r;
        (*liniarizare)[i][e + 11].g = culoare.g;
        (*liniarizare)[i][e + 11].b = culoare.b;
    }


}

int comparareDetectii(const void *detectie1, const void *detectie2) {
    if (((struct fi *) detectie1)->correlation > ((struct fi *) detectie2)->correlation) return -1;
    if (((struct fi *) detectie1)->correlation == ((struct fi *) detectie2)->correlation) return 0;
    return 1;

}



//todo functii CRIPTARE DECRIPTARE

unsigned int xorshift32(unsigned int seed) {

    unsigned int r = seed;

    r = r ^ r << 13;
    r = r ^ r >> 17;
    r = r ^ r << 5;

    return r;

}


void CitireLiniarizare(char *bmpIn, unsigned char **header, struct px **liniarizare, unsigned int *latime_img,
                       unsigned int *inaltime_img) {

    FILE *fin;


    fin = fopen(bmpIn, "rb");

    if (fin == NULL) {
        printf("nu am gasit imaginea sursa din care citesc");
        return;
    }

    //citesc dimensiunile de la octetul 18 respectiv 22
    fseek(fin, 18, SEEK_SET);
    fread(latime_img, sizeof(unsigned int), 1, fin);
    fread(inaltime_img, sizeof(unsigned int), 1, fin);
    printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n", (*latime_img), (*inaltime_img));
    rewind(fin);

    int k;

    *header = (unsigned char *) malloc(54 * sizeof(unsigned char));
    *liniarizare = (struct px *) malloc(((3 * (*inaltime_img) * (*latime_img))) * sizeof(unsigned char));


    fread(*header, 1, 54, fin);


    fseek(fin, 54, SEEK_SET);
    int padding= (4- (*latime_img*3)%4)%4;

    int i, j;
    for (i = 0; i < *inaltime_img; i++) {
        int var_ajutatoare=(*inaltime_img - i - 1);
        int ct = *latime_img * var_ajutatoare;
        for (j = 0; j < *latime_img; j++) {
            fread((*liniarizare + ct), 1, 3, fin);
            ct++;
        }
        //todo nu uita de padding
        fseek(fin,padding *sizeof(char), SEEK_CUR );
    }


    fclose(fin);

}

void writeBmpOut1(struct px *liniarizare, char *bmpOut1, char *bmpIn, unsigned char *header, unsigned int latime_img,
                  unsigned int inaltime_img) {

    int i;
    FILE *fout = fopen(bmpOut1, "wb");
    rewind(fout);
    //scriu headerul in fisierul binar
    fwrite(header, 1, 54, fout);

    int octNecesari= (4- (latime_img*3)%4)%4;
    int ct, j;
    for (i = 0; i < inaltime_img; i++) {
        int var_ajutatoare=(inaltime_img - i - 1);
        ct = latime_img * var_ajutatoare;
        for (j = 0; j < latime_img; ++j) {
            fwrite((liniarizare + ct), sizeof(char), 3, fout);
            ++ct;

        }
        //todo nu uita de padding
        fwrite(liniarizare+ct-1,sizeof(char),octNecesari,fout);
    }
    fclose(fout);
}

void criptare(char *bmpIn, char *bmpOut, char *secretKey, struct px **liniarizare, unsigned int latime_img,
              unsigned int inaltime_img, unsigned int *R0, unsigned int SV) {

    unsigned int i;


    /* //dechid fisierul sursa
     FILE *fout = fopen(bmpOut, "rb+");
     if (fout == NULL) {
         printf("nu am gasit imaginea sursa din care citesc");
         return;
     }*/



    //generez W*H-1 nr aleatoare cu xorshift
    unsigned int dimR = (inaltime_img * latime_img);




    //initializez permutarea
    unsigned int *P, aux;
    P = (unsigned int *) malloc((dimR) * sizeof(unsigned int));
    for (i = 0; i < dimR; i++)
        P[i] = i;


    unsigned int r,j = 0;

    r=xorshift32(*R0);
    for (i = dimR-1; i >= 1; i--) {
        j =r%(i+1);

        aux = P[i];
        P[i] = P[j];
        P[j] = aux;

        r=xorshift32(r);
    }

    //efectuez permutarea efectiva pe poza
    struct px *auxx;
    auxx = (struct px *) malloc((dimR) * sizeof(struct px));

    for (i = 0; i < dimR; i++)
        auxx[P[i]] = *(*liniarizare+i);
    for (i = 0; i < dimR; i++)
        *(*liniarizare+i)=auxx[i];


    *R0=r;

    //xorez elementele


}

void xorare1(struct px **liniarizare,unsigned int SV,unsigned int R0,unsigned int latime_img,unsigned int inaltime_img){
    int i=0;

    while(i<(latime_img*inaltime_img)) {

        if(i==0) {
            //xorez primul element
            unsigned int z = 0;
            unsigned int copie1=SV;
            z = copie1 & ((1 << 8) - 1);
            copie1 = copie1 >> 8;
            (*liniarizare)[i].b = (*liniarizare)[i].b ^ z;

            z = copie1 & ((1 << 8) - 1);
            copie1 = copie1 >> 8;
            (*liniarizare)[i].g = (*liniarizare)[i].g ^ z;

            z = copie1 & ((1 << 8) - 1);
            copie1 = copie1 >> 8;
            (*liniarizare)[i].r = (*liniarizare)[i].r ^ z;

            z = 0;
            unsigned int copie2=R0;
            z = copie2 & ((1 << 8) - 1);
            copie2 = copie2 >> 8;
            (*liniarizare)[i].b = (*liniarizare)[i].b ^ z;

            z = copie2 & ((1 << 8) - 1);
            copie2 = copie2 >> 8;
            (*liniarizare)[i].g = (*liniarizare)[i].g ^ z;

            z = copie2 & ((1 << 8) - 1);
            copie2 = copie2 >> 8;
            (*liniarizare)[i].r = (*liniarizare)[i].r ^ z;

            R0 = xorshift32(R0);
            i++;
        }
        else{
            unsigned int z = 0;
            unsigned int copie1=SV;
            unsigned int copie2=R0;
            (*liniarizare)[i].g = (*liniarizare)[i].g ^ (*liniarizare)[i - 1].g;
            (*liniarizare)[i].r = (*liniarizare)[i].r ^ (*liniarizare)[i - 1].r;
            (*liniarizare)[i].b = (*liniarizare)[i].b ^ (*liniarizare)[i - 1].b;

            z = 0;
            z = copie2 & ((1 << 8) - 1);
            copie2 = copie2 >> 8;
            (*liniarizare)[i].b = (*liniarizare)[i].b ^ z;

            z = copie2 & ((1 << 8) - 1);
            copie2 = copie2 >> 8;
            (*liniarizare)[i].g = (*liniarizare)[i].g ^ z;

            z = copie2 & ((1 << 8) - 1);
            copie2 = copie2 >> 8;
            (*liniarizare)[i].r = (*liniarizare)[i].r ^ z;

            copie2 = xorshift32(copie2);
            R0=xorshift32(R0);
            i++;
        }

    }

}

void chiPatrat(struct px *liniarizare, unsigned int latime_img, unsigned int inaltime_img) {

    FILE *fout = fopen("chiTestValues.txt", "w");
    double frecv2, chi2 = 0;
    unsigned int *F, i;
    frecv2 = (latime_img * inaltime_img) / 256;
    F = (unsigned int *) calloc(256, sizeof(unsigned int));


    //afisez pt canalul B
    for (i = 0; i < latime_img * inaltime_img; i++)
        F[((liniarizare)[i]).b]++;

    for (i = 0; i < 255; i++) {
        chi2 = chi2 + ((F[i] - frecv2) * (F[i] - frecv2)) / frecv2;

    }
    fprintf(fout, "\nValorile testului chi patrat  pentru enc_peppers1.bmp sunt:");
    fprintf(fout, "\n B: %.2lf", chi2);


    //pt canalul G
    for (i = 0; i <= 255; i++)F[i] = 0;
    chi2 = 0;
    for (i = 0; i < latime_img * inaltime_img; i++)
        F[((liniarizare)[i]).g]++;

    for (i = 0; i < 255; i++) {
        chi2 = chi2 + ((F[i] - frecv2) * (F[i] - frecv2)) / frecv2;

    }
    fprintf(fout, "\n\n G: %.2lf", chi2);


    //pentru canalul R
    for (i = 0; i <= 255; i++)F[i] = 0;
    chi2 = 0;
    for (i = 0; i < latime_img * inaltime_img; i++)
        F[((liniarizare)[i]).r]++;

    for (i = 0; i < 255; i++) {
        chi2 = chi2 + ((F[i] - frecv2) * (F[i] - frecv2)) / frecv2;

    }
    fprintf(fout, "\n\n R: %.2lf", chi2);
}

void decriptare(char *bmpIn, char *bmpOut,struct px ** liniarizare,unsigned int latime_img,unsigned int inaltime_img, unsigned int *R0,unsigned int SV){
    unsigned int i,*P,*PInv;
    unsigned int dimR = (inaltime_img * latime_img);
    P=(unsigned int* )malloc (dimR* sizeof(unsigned int));
    for (i = 0; i < dimR; i++)
        P[i] = i;


    unsigned int r,j = 0,aux;

    r=xorshift32(*R0);
    for (i = dimR-1; i >= 1; i--) {
        j =r%(i+1);

        aux = P[i];
        P[i] = P[j];
        P[j] = aux;

        r=xorshift32(r);
    }


    //creez permutarea inversa
    PInv=(unsigned int* )malloc (dimR* sizeof(unsigned int));
    for (i = 0; i < dimR; i++)
        PInv[P[i]]=i;
    xorare1(&(*liniarizare),SV,*R0,latime_img,inaltime_img);


    //permut dupa inversa
    struct px *auxx;
    auxx = (struct px *) malloc((dimR) * sizeof(struct px));

    for (i = 0; i < dimR; i++)
        auxx[PInv[i]] = *(*liniarizare+i);
    for (i = 0; i < dimR; i++)
        *(*liniarizare+i)=auxx[i];


    *R0=r;

}

int main() {


    char img[100],s0[50],s1[50],s2[50],s3[50],s4[50],s5[50],s6[50],s7[50],s8[50],s9[50];
    FILE *fin=fopen("sabloane.txt","r");
    fscanf(fin,"%s",img);
    fscanf(fin,"%s",s0);
    fscanf(fin,"%s",s1);
    fscanf(fin,"%s",s2);
    fscanf(fin,"%s",s3);
    fscanf(fin,"%s",s4);
    fscanf(fin,"%s",s5);
    fscanf(fin,"%s",s6);
    fscanf(fin,"%s",s7);
    fscanf(fin,"%s",s8);
    fscanf(fin,"%s",s9);



    char imgOut[] = "testOut.bmp";



    struct px **main_img, **cifra0, *culoare, **cifra1, **cifra2, **cifra3, **cifra4, **cifra5, **cifra6, **cifra7, **cifra8, **cifra9;
    unsigned int H_main, W_main, H_sablon, W_sablon, n = 165, numarDetectii = 0, contor_culoare,q=0,contorSuprapuneri;
    unsigned char *header2, *header3;
    double ps = 0.5;
    int contor1;
    struct fi *D, fiCurent,*Dsuprapuneri;
    D = (struct fi *) malloc(5000 * sizeof(struct fi));
    Dsuprapuneri = (struct fi *) malloc(1000 * sizeof(struct fi));


    incarcareImg2(img, &H_main, &W_main, &main_img, &header2);
    grayscale(img, H_main, W_main, &main_img);
    writeBmpOut2(main_img, imgOut, img, header2, W_main, H_main);

    incarcareImg2(s0, &H_sablon, &W_sablon, &cifra0, &header3);
    grayscale(s0, H_sablon, W_sablon, &cifra0);

    incarcareImg2(s1, &H_sablon, &W_sablon, &cifra1, &header3);
    grayscale(s1, H_sablon, W_sablon, &cifra1);

    incarcareImg2(s2, &H_sablon, &W_sablon, &cifra2, &header3);
    grayscale(s2, H_sablon, W_sablon, &cifra2);

    incarcareImg2(s3, &H_sablon, &W_sablon, &cifra3, &header3);
    grayscale(s3, H_sablon, W_sablon, &cifra3);

    incarcareImg2(s4, &H_sablon, &W_sablon, &cifra4, &header3);
    grayscale(s4, H_sablon, W_sablon, &cifra4);

    incarcareImg2(s5, &H_sablon, &W_sablon, &cifra5, &header3);
    grayscale(s5, H_sablon, W_sablon, &cifra5);

    incarcareImg2(s6, &H_sablon, &W_sablon, &cifra6, &header3);
    grayscale(s6, H_sablon, W_sablon, &cifra6);

    incarcareImg2(s7, &H_sablon, &W_sablon, &cifra7, &header3);
    grayscale(s7, H_sablon, W_sablon, &cifra7);

    incarcareImg2(s8, &H_sablon, &W_sablon, &cifra8, &header3);
    grayscale(s8, H_sablon, W_sablon, &cifra8);

    incarcareImg2(s9, &H_sablon, &W_sablon, &cifra9, &header3);
    grayscale(s9, H_sablon, W_sablon, &cifra9);




    printf("\n-------------Desenare  dreptunghuri -----------\n");


    initializareCulori(&culoare);

    corelatie(main_img, cifra0, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra1, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra2, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra3, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra4, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra5, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra6, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra7, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra8, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);
    corelatie(main_img, cifra9, ps, &D, &fiCurent, H_main, W_main, H_sablon, W_sablon, &numarDetectii,&q);



    qsort(D, numarDetectii + 2, sizeof(struct fi), comparareDetectii);



    //Afisari pe ecran
    printf("\n---Printam corelatiile detectiilor sortate ---\n");

    for (contor1 = 1; contor1 <= numarDetectii; contor1++)
        printf("%.2f; ", (D)[contor1].correlation);


    printf("-------Numar total de detectii----------: %d\n", numarDetectii);

    printf("\n---Desenam dreptunghuri ---\n");
    for (contor1 = 1; contor1 < numarDetectii; contor1++) {
        printf("Desenam dreptunghiul cu numarul %d...\n", contor1);
        desenareDreptunghiuri(&main_img, D[contor1], culoare[D[contor1].culoare], numarDetectii);
    }



    //eliminareMaxime(D,&main_img,culoare,numarDetectii,H_sablon,W_sablon,&contorSuprapuneri,&Dsuprapuneri);
    writeBmpOut2(main_img, imgOut, img, header2, W_main, H_main);


    char bmpIn[] = "peppers.bmp";
    char secretKey[] = "secret_key.txt";
    char bmpOut1[] = "enc_peppers1.bmp";
    unsigned char *header;
    struct px *liniarizare;//creez vectorul pentru liniarizarea matricei
    int i;
    unsigned int latime_img, inaltime_img;
    unsigned int R0;
    unsigned int SV;

    FILE *fsecret = fopen(secretKey, "r");
    if (fsecret == NULL) {
        printf("nu am gasit secret_key.txt din care sa citesc");
        return 0;
    }
    //citesc valorile din cheia secreta

    fscanf(fsecret, "%u", &R0);

    fscanf(fsecret, "%u", &SV);


    CitireLiniarizare(bmpIn, &header, &liniarizare, &latime_img, &inaltime_img);
    criptare(bmpIn, bmpOut1, secretKey, &liniarizare, latime_img, inaltime_img,&R0,SV);
    xorare1(&liniarizare,SV,R0,latime_img,inaltime_img);
    //decriptare(bmpIn, bmpOut1,&liniarizare,latime_img,inaltime_img,&R0,SV);
    writeBmpOut1(liniarizare, bmpOut1, bmpIn, header, latime_img, inaltime_img);
    chiPatrat(liniarizare, latime_img, inaltime_img);

    fclose(fsecret);
    free(liniarizare);
    free(header);
    free(header2);
    free(header3);
    free(D);
    free(Dsuprapuneri);
    free(culoare);
    fclose(fin);
    return 0;
}
