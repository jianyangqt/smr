#ifndef MY_FUNCTIONS
#define MY_FUNCTIONS

#define KEY_LEN 128
struct Key_linker{
    char key[KEY_LEN];
    void * value;
    struct Key_linker * next;

};


unsigned short  add_hash_key_value(struct Key_linker *[], const char *, struct Key_linker *);
struct Key_linker * get_value_of_key(struct Key_linker *[], const char *);
int check_key_exist(struct Key_linker *[], const char *);

extern unsigned short hash_func(const char *);
#endif
