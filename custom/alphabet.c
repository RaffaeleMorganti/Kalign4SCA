/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "tldevel.h"
#include <stdint.h>

#define ALPHABET_IMPORT
#include "alphabet.h"


int create_default_protein(struct alphabet* a);
int create_default_DNA(struct alphabet* a);
int create_reduced_protein(struct alphabet *a);

int clean_and_set_to_extern(struct alphabet* a);

static int merge_multiple(struct alphabet*a,char* p,int n);
static int merge_codes(struct alphabet*a,const int X, const int Y);

#ifdef UTEST_ALPHABET
int print_alphabet(struct alphabet* a);


int main(void)
{
        struct alphabet* a = NULL;

        RUNP(a = create_alphabet(ALPHA_defPROTEIN));

        print_alphabet(a);
        MFREE(a);
        a = NULL;
        RUNP(a = create_alphabet(ALPHA_redPROTEIN));

        print_alphabet(a);
        MFREE(a);

        RUNP(a = create_alphabet(ALPHA_defDNA));

        print_alphabet(a);
        MFREE(a);

        RUNP(a = create_alphabet(ALPHA_redPROTEIN2));

        print_alphabet(a);
        MFREE(a);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}



int print_alphabet(struct alphabet* a)
{
        fprintf(stdout,"LEN: %d\n",a->L);
        int i;
        for(i = 64;i < 96;i++){

                fprintf(stdout,"%c\t%d\n",  (char)i, a->to_internal[i]);
        }
        return OK;
}
#endif

struct alphabet* create_alphabet(int type)
{
        struct alphabet* a = NULL;
        int i;
        MMALLOC(a, sizeof(struct alphabet));

        for(i = 0; i < 128;i++){
                a->to_internal[i] = -1;
        }

        for(i = 0; i < 32;i++){
                a->to_external[i] = -1;
        }

        switch (type) {
        case ALPHA_defPROTEIN :
        case ALPHA_ambigiousPROTEIN :{
                create_default_protein(a);
                break;
        }
        case ALPHA_defDNA : {
                create_default_DNA(a);
                break;
        }
        case ALPHA_redPROTEIN :
        case ALPHA_redPROTEIN2 : {
                create_reduced_protein(a);
                break;
        }
        default:
                break;
        }
        RUN(clean_and_set_to_extern(a));
        return a;
ERROR:
        if(a){
                MFREE(a);
        }
        return NULL;
}

int switch_alphabet(struct alphabet* a, int type)
{
        int i;
        for(i = 0; i < 128;i++){
                a->to_internal[i] = -1;
        }
        for(i = 0; i < 32;i++){
                a->to_external[i] = -1;
        }

        switch (type) {
        case ALPHA_defPROTEIN : {
                create_default_protein(a);
                break;
        }
        case ALPHA_redPROTEIN : {
                create_reduced_protein(a);
                break;
        }
        default:
                break;
        }

        RUN(clean_and_set_to_extern(a));
        return OK;
ERROR:
        return FAIL;
}


int create_default_protein(struct alphabet* a)
{
        char aacode[18] = "ABCDEFGHIJKLMNOPQR";

        int code;
        int i;
        code = 0;
        for(i = 0; i < 18;i++){
                //fprintf(stdout,"%c %d CODE: %d\n", aacode[i], (int) aacode[i], code);
                a->to_internal[(int) aacode[i]] = code;

                code++;
        }
        
        return OK;
}


int create_default_DNA(struct alphabet* a)
{
        char dnacode[6] = "UVWXYZ";

        int code;
        int i;
        code = 0;
        for(i = 0; i < 6;i++){
                //fprintf(stdout,"%c %d CODE: %d\n", aacode[i], (int) aacode[i], code);
                a->to_internal[(int) dnacode[i]] = code;

                code++;
        }
        
        return OK;
}


int create_reduced_protein(struct alphabet* a)
{
        char aacode[18] = "ABCDEFGHIJKLMNOPQR";

        int code;
        int i;
        code = 0;
        for(i = 0; i < 18;i++){
                a->to_internal[(int) aacode[i]] = code;
                code++;
        }

        /* reduced codes */

        merge_multiple(a,"ABC",3);
        merge_multiple(a,"DEF",3);
        merge_multiple(a,"GHI",3);
        merge_multiple(a,"JKL",3);
        merge_multiple(a,"MNO",3);
        merge_multiple(a,"PQR",3);
        
        return OK;

}


int merge_multiple(struct alphabet*a,char* p,int n)
{

 
        // Declaring pointer to the
        // argument list
        int min = INT32_MAX;

        for(int i = 0; i < n;i++){
                min = MACRO_MIN(min,a->to_internal[(int)p[i]]);
        }

        for(int i = 0; i < n;i++){
                a->to_internal[(int)p[i]] = min;
        }
        return OK;
}


int merge_codes(struct alphabet*a,const int X, const int Y)
{
        int min;

        min = MACRO_MIN(a->to_internal[X],a->to_internal[Y]);

        ASSERT(min != -1, "code not set!");

        a->to_internal[X] = min;
        a->to_internal[Y] = min;
        return OK;
ERROR:
        return FAIL;
}



int clean_and_set_to_extern(struct alphabet* a)
{
        int i;
        int code = 0;
        int8_t trans[32];
        for(i = 0; i < 32;i++){
                trans[i] = -1;

        }

        for(i = 64; i < 96;i++){
                if(a->to_internal[i] != -1){
                        trans[a->to_internal[i]] = 1;
                }
        }
        code = 0;
        for(i = 0; i < 32;i++){
                if(trans[i] == 1){
                        trans[i] = code;
                        code++;
                }
        }
        a->L = code;
        for(i = 64; i < 96;i++){
                if(a->to_internal[i] != -1){
                        a->to_internal[i] = trans[a->to_internal[i]];//a->to_internal[i]];
                        a->to_internal[i+32] = a->to_internal[i];

                }

        }

        for(i = 64;i < 96;i++){
                if(a->to_internal[i] != -1){
                        a->to_external[a->to_internal[i]] = i;
                }
        }
        return OK;
}
