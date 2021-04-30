#include <string.h>
#include <inttypes.h>
#include <stdlib.h>

#define KEY_LEN 128
struct Key_linker{
    char key[KEY_LEN];
    void * value;
    struct Key_linker * next;

};



//Hash function using short.
unsigned short
hash_func(const char * short_str)
{
    unsigned short dt_out = 0;
    unsigned char j = 0;
    unsigned int k = 0;
    unsigned short peirod_amplifer = 170;
    for (int i = 0; i < strlen(short_str); i++){
        if ((k % 3) == 0){
            if (j == 0){
                j = 1;
                dt_out ^= (unsigned short)((short_str[i] | peirod_amplifer) << 8);

            } else{
                j = 0;
                dt_out ^= short_str[i] | peirod_amplifer;
            }

        } else if ((k % 3) == 1){
            if (j == 0){
                j = 1;
                dt_out ^= (unsigned short)((short_str[i] & peirod_amplifer) << 8);

            } else{
                j = 0;
                dt_out ^= short_str[i] & peirod_amplifer;
            }

        } else if ((k % 3) == 2){
            if (j == 0){
                j = 1;
                dt_out ^= (unsigned short)((short_str[i] ^ peirod_amplifer) << 8);

            } else{
                j = 0;
                dt_out ^= short_str[i] ^ peirod_amplifer;
            }

        }
        
        k++;
    }

    return dt_out;
}



unsigned short
add_hash_key_value(struct Key_linker ** hash_table, const char * key, struct Key_linker * key_linker)
{
    unsigned short hash_value = hash_func(key);
    struct Key_linker * tmp;

    if (!hash_table[hash_value]){
        hash_table[hash_value] = key_linker;

    } else{
        tmp = hash_table[hash_value];
        while(tmp -> next){

            tmp = tmp -> next;
        }

        tmp -> next = key_linker;
    }
    
    return hash_value;
}


int
check_key_exist(struct Key_linker ** hash_table, const char * key){
    unsigned short hash_value = hash_func(key);
    struct Key_linker * linker;
    linker = hash_table[hash_value];
    if(!linker){
        return 0;

    } else{
        while(linker){
            if (strcmp(key, linker -> key) == 0){
                return 1;
            }
            linker = linker -> next; 
        }

    }
    
    return 0;
}



struct key_linker *
get_value_of_key(struct key_linker ** hash_table, const char * key)
{
    struct key_linker * linker;
    unsigned short hash_value = hash_func(key);
    linker = hash_table[hash_value];
    return linker;
}






