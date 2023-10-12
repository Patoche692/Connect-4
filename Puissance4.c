#define _CRT_SECURE_NO_WARNINGS
#include <stdint.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define P4_COLONNES (7)
#define P4_LIGNES (6)
#define HEIGHT (6)
#define J1_JETON ('O')
#define J2_JETON ('X')

#define ACT_ERR (0)
#define ACT_JOUER (1)
#define ACT_NOUVELLE_SAISIE (2)
#define ACT_QUITTER (3)

#define STATUT_OK (0)
#define STATUT_GAGNE (1)
#define STATUT_EGALITE (2)
#define PROFONDEUR (5)
#define WIDTH (7)  // width of the board
#define HEIGHT (6) // height of the board
#define MIN_SCORE (-(WIDTH*HEIGHT)/2 + 3)
#define MAX_SCORE ((WIDTH*HEIGHT+1)/2 - 3)

struct position
{
    int colonne;
    int ligne;
};

struct Entry {
  long long unsigned int key : 56; // use 56-bit keys
  char val;      // use 8-bit values
};  

struct MoveScore
{
  long long unsigned int move;
  int score;
};

static void affiche_grille(int jeton);
static int demande_action(int *);
static int demande_nb_joueur(void);
double nb_aleatoire(void);
int nb_aleatoire_entre(int, int);
static int position_valide(struct position *);
static int statut_jeu();
int IA_jouer();
int gagnant();
int nb_series(long long unsigned int pos);
int alphabeta(int profondeur, int alpha, int beta, int lastCol);
void printb(long long int n);
long long unsigned int alignment(long long unsigned int pos, int n);
int canPlay(int col);
static long long unsigned int top_mask(int col);
void playCol(int col);
void play(long long unsigned int move);
static long long unsigned int bottom_mask_col(int col);
static long long unsigned int column_mask(int col);
void change_current_pos();
unsigned int entry_index(long long unsigned int key);
void put(long long unsigned int key, char val);
char get(long long unsigned int key);
void reset();
int negamaxStart(int profondeur, int alpha, int beta);
int negamaxEnd(int alpha, int beta);
int canWinNext();
long long unsigned int possible();
int isWinningMove(int col);
long long unsigned int winning_position();
long long unsigned int opponent_winning_position();
static long long unsigned int compute_winning_position();
long long unsigned int possibleNonLoosingMoves();
void init();
static long long unsigned int bottom(int width, int height);
static long long unsigned int compute_winning_position(long long unsigned int position);
void add(long long unsigned int move, int score, struct MoveScore* moveScore, int* size);
long long unsigned int getNext(struct MoveScore* moveScore, int* size);
int move_score(long long unsigned int move);
int solve(int profondeur);
static long long unsigned int compute_alignment_three_position(long long unsigned int position);
static unsigned int popcount(long long unsigned int m);
int move_score_end(long long unsigned int move);

int lastSeries1 = 0, lastSeries2 = 0;
int nbJetons = 0, moves = 0;
long long unsigned int currentPosition = 0b0000000000000000000000000000000000000000000000000000000000000000,
mask = 0b0000000000000000000000000000000000000000000000000000000000000000, 
test2 = 0b0000000000000000000000000000000000000000000000000000000000000000;
struct Entry entries[8388593];
long long unsigned int bottom_mask;
long long unsigned int board_mask;
int nodes = 0;

void init(){
  bottom_mask = bottom(WIDTH, HEIGHT);
  board_mask = bottom_mask * ((1LL << HEIGHT)-1);
}

static void affiche_grille(int jeton)
{
    /*
     * Affiche la grille pour le ou les joueurs.
     */

    int col;
    int lgn;
    long long unsigned int tmp;
    putchar('\n');

    for (col = 1; col <= P4_COLONNES; ++col)
        printf("  %d ", col);

    putchar('\n');
    putchar('+');

    for (col = 1; col <= P4_COLONNES; ++col)
        printf("---+");

    putchar('\n');

    for (lgn = 0; lgn < P4_LIGNES; ++lgn)
    {
        putchar('|');
        tmp = (1LL << 5) >> lgn;

        for (col = 0; col < P4_COLONNES; ++col){
            if (mask & tmp){
              if(jeton == J1_JETON){
                if((currentPosition & tmp)){
                  printf(" O |");
                }
                else{
                  printf(" X |");
                }
              }
              else{
                if((currentPosition & tmp)){
                  printf(" X |");
                }
                else{
                  printf(" O |");
                }
              }
            }
            else
                printf(" %c |", ' ');


            tmp <<= 7;
        }

        putchar('\n');
        putchar('+');

        for (col = 1; col <= P4_COLONNES; ++col)
            printf("---+");

        putchar('\n');
    }

    for (col = 1; col <= P4_COLONNES; ++col)
        printf("  %d ", col);

    putchar('\n');
}

static int demande_action(int *coup)
{
    /*
     * Demande l'action à effectuer au joueur courant.
     * S'il entre un chiffre, c'est qu'il souhaite jouer.
     * S'il entre la lettre « Q » ou « q », c'est qu'il souhaite quitter.
     * S'il entre autre chose, une nouvelle saisie sera demandée.
     */

    char c;
    int ret = ACT_ERR;

    if (scanf("%d", coup) != 1)
    {
        if (scanf("%c", &c) != 1)
        {
            fprintf(stderr, "Erreur lors de la saisie\n");
            return ret;
        }

        switch (c)
        {
        case 'Q':
        case 'q':
            ret = ACT_QUITTER;
            break;
        default:
            ret = ACT_NOUVELLE_SAISIE;
            break;
        }
    }
    else
        ret = ACT_JOUER;

    return ret;
}


static int demande_nb_joueur(void)
{
    /*
     * Demande et récupère le nombre de joueurs.
     */

    int njoueur = 0;

    while (1)
    {
        printf("Combien de joueurs prennent part à cette partie ? ");

        if (scanf("%d", &njoueur) != 1 && ferror(stdin))
        {
                fprintf(stderr, "Erreur lors de la saisie\n");
                return 0;
        }
        else if (njoueur != 1 && njoueur != 2)
            fprintf(stderr, "Plait-il ?\n");
        else
            break;
    }

    return njoueur;
}


double nb_aleatoire(void)
{
    /*
     * Retourne un nombre pseudo-aléatoire compris entre zéro inclus et un exclus.
     */

    return rand() / ((double)RAND_MAX + 1.);
}


int nb_aleatoire_entre(int min, int max)
{
    /*
     * Retourne un nombre pseudo-aléatoire entre `min` et `max` inclus.
     */

    return nb_aleatoire() * (max - min + 1) + min;
}


static int statut_jeu()
{
    /*
     * Détermine s'il y a lieu de continuer le jeu ou s'il doit être
     * arrêté parce qu'un joueur a gagné ou que la grille est complète.
     */

    if (nbJetons == 42)
        return STATUT_EGALITE;
    if (alignment(currentPosition, 4))
        return STATUT_GAGNE;

    return STATUT_OK;
}


void printb(long long int n)
{
	long long int bit = 0 ;
	long long int mask = 1 ;
	for (int i = 0 ; i < 64 ; ++i)
	{
		bit = (n & mask) >> i ;
		printf("%lld", bit) ;
		mask <<= 1 ;
	}
}

int main(void)
{
    int statut;
    char jeton = J1_JETON;
    int njoueur;
    srand(time(NULL));
    init();
    affiche_grille(jeton);
    njoueur = demande_nb_joueur();

    if (!njoueur)
        return EXIT_FAILURE;

    while (1)
    {
        struct position pos;
        int action;
        int coup;

        if (njoueur == 1 && jeton == J2_JETON)
        {

            coup = IA_jouer();
            printf("Joueur 2 : %d\n", coup + 1);
        }
        else
        {
            printf("Joueur %d : ", (jeton == J1_JETON) ? 1 : 2);
            action = demande_action(&coup);

            if (action == ACT_ERR)
                return EXIT_FAILURE;
            else if (action == ACT_QUITTER)
                return 0;
            else if (action == ACT_NOUVELLE_SAISIE || !canPlay(coup - 1))
            {
                fprintf(stderr, "Vous ne pouvez pas jouer à cet endroit\n");
                continue;
            }
            coup--;
        }

        playCol(coup);
        printf("\n");
        nbJetons++;
        affiche_grille(jeton);
        statut = statut_jeu();

        if (statut != STATUT_OK)
            break;

        jeton = (jeton == J1_JETON) ? J2_JETON : J1_JETON;
    }

    if (statut == STATUT_GAGNE)
        printf("Le joueur %d a gagné\n", (jeton == J1_JETON) ? 1 : 2);
    else if (statut == STATUT_EGALITE)
        printf("Égalité\n");

    return 0;
}

int IA_jouer(){
  time_t lastTime = time(NULL);
  int profondeur = 1, best = 0;
  if(nbJetons < 15)
    while(time(NULL) - lastTime < 2 && profondeur < 42 - nbJetons){
      nodes = 0;
      best = solve(profondeur+=2);
      reset();
    }
  else{
    nodes = 0;
    best = solve(0);
    profondeur = 42 - nbJetons;
  }
  printf("Profondeur : %d\n", profondeur);
  printf("Nombre de noeuds : %d\n", nodes);
  return best;
}

int solve(int profondeur)
{
	int max = -129 + nbJetons;
	int alpha = -128 + nbJetons;
	int beta = 128 - nbJetons;
	int tmp;
	int bestMove[7] = {0}, nbBestMove = -1;
  long long unsigned int lastPosition = currentPosition, lastMask = mask;

  for(int i = 0; i < 4; i=-i + (-i >= 0 ? 1 : 0)){
    if(canPlay(3 + i)){
      nbJetons++;
      playCol(3 + i);
      if(gagnant() || nbJetons == 42){
        nbJetons--;
        currentPosition = lastPosition;
        mask = lastMask;
        return 3 + i;
      }
        nbJetons--;
        currentPosition = lastPosition;
        mask = lastMask;
    }
  }

  long long unsigned int next = possibleNonLoosingMoves();

  if(!profondeur){
    int bestScore[7] = {0};
    for(int i = 0; i < 4; i=-i + (-i >= 0 ? 1 : 0))
    {
      if(next & column_mask(3 + i)){
        nodes++;
        playCol(3 + i);
        nbJetons++;

        tmp = -negamaxEnd(-beta, -alpha);
        if(tmp >= max)
        {
          if(tmp > max){
            max = tmp;
            nbBestMove = -1;
          }
          nbBestMove++;
          bestMove[nbBestMove] = 3 + i;
          bestScore[nbBestMove] = tmp;
        }
        currentPosition = lastPosition;
        mask = lastMask;
        nbJetons--;
      }
      else if(!next && canPlay(3 + i)){
        return 3 + i;
      }

    }
    if(bestScore[0] > 50)
      printf("Vous perdez dans %d coups\n", 128 - nbJetons - bestScore[0]);
    else if(bestScore[0] < 50){
      printf("Vous pouvez gagner dans %d coups\n", 128 - nbJetons + bestScore[0]);
    }
    return bestMove[nb_aleatoire_entre(0, nbBestMove)];
  }

	for(int i = 0; i < 4; i=-i + (-i >= 0 ? 1 : 0))
	{
    if(next & column_mask(3 + i)){
      nodes++;
  		playCol(3 + i);
  		nbJetons++;

  		tmp = -negamaxStart(profondeur-1, -beta, -alpha);
  		if(tmp >= max)
  		{
  			if(tmp > max){
  				max = tmp;
  				nbBestMove = -1;
  			}
  			nbBestMove++;
  			bestMove[nbBestMove] = 3 + i;
  			
  		}
      currentPosition = lastPosition;
      mask = lastMask;
  		nbJetons--;
  	}
    else if(!next && canPlay(3 + i)){
      return 3 + i;
    }

  }

	return bestMove[nb_aleatoire_entre(0, nbBestMove)];
}

int gagnant()
{
  if(alignment(currentPosition, 4)){
    return 1;
  }
  return 0;
}

int negamaxStart(int profondeur, int alpha, int beta) {

  int score = 0, val = 0, max = 128 - nbJetons;

  long long unsigned int possible = possibleNonLoosingMoves();

  if(!possible){
    return -128 + nbJetons;
  }

  if(!profondeur)
  {
    return /*popcount(compute_winning_position(currentPosition ^ mask)) - popcount(compute_winning_position(currentPosition));*/nb_series(currentPosition ^ mask) - nb_series(currentPosition);
  }

  if(val = get(currentPosition + mask)){
    max = val;
  }
  if(beta > max){
    beta = max;
    if(alpha >= beta) return beta;
  }

  struct MoveScore moveScore[WIDTH];
  int size = 0;
  long long unsigned int move;
  for(int i = 0; i < 4; i=-i + (-i >= 0 ? 1 : 0))
    if((move = possible & column_mask(3 + i)))
      add(move, move_score(move), moveScore, &size);

  long long unsigned int next = 0, lastPosition = currentPosition, lastMask = mask;
  /*for(int i = 0; i < 4; i=-i + (-i >= 0 ? 1 : 0)) 
    if(possible & column_mask(3 + i)){*/
    while((next = getNext(moveScore, &size))){
      nodes++;
      play(next);
      nbJetons++;
      score = -negamaxStart(profondeur - 1, -beta, -alpha);
      currentPosition = lastPosition;
      mask = lastMask;
      nbJetons--;

      if(score >= beta) return score;
      if(score > alpha) alpha = score;
    }

  put(currentPosition + mask, alpha);
  return alpha;
}

int negamaxEnd(int alpha, int beta){

  int score = 0, val = 0, max = 128 - nbJetons;

  long long unsigned int possible = possibleNonLoosingMoves();

  if(!possible){
    return -128 + nbJetons;
  }

  if(val = get(currentPosition + mask)){
    max = val;
  }
  if(beta > max){
    beta = max;
    if(alpha >= beta) return beta;
  }

  struct MoveScore moveScore[WIDTH];
  int size = 0;
  long long unsigned int move;
  for(int i = 0; i < 4; i=-i + (-i >= 0 ? 1 : 0))
    if((move = possible & column_mask(3 + i)))
      add(move, move_score_end(move), moveScore, &size);

  long long unsigned int next = 0, lastPosition = currentPosition, lastMask = mask;
  while((next = getNext(moveScore, &size))){
    nodes++;
    play(next);
    nbJetons++;
    score = -negamaxEnd(-beta, -alpha);
    currentPosition = lastPosition;
    mask = lastMask;
    nbJetons--;

    if(score >= beta) return score;
    if(score > alpha) alpha = score;
  }

  put(currentPosition + mask, alpha);
  return alpha;
}

int nb_series(long long unsigned int pos)
//Compte le nombre de séries de n pions alignés de chacun desjoueurs
{
  int tmp = 0;
  long long unsigned int d1, d2, v, h = pos & (pos >> (HEIGHT+1));
  if(h = (h & (h >> ((HEIGHT+1))))) tmp = 1;

  // diagonal 1
  d1 = pos & (pos >> HEIGHT);
  if(d1 = (d1 & (d1 >> (HEIGHT)))) tmp = 1;

  // diagonal 2 
  d2 = pos & (pos >> (HEIGHT+2));
  if(d2 = (d2 & (d2 >> ((HEIGHT+2))))) tmp = 1;

  // vertical;
  v = pos & (pos >> 1);
  if(v = (v & (v >> 1))) tmp = 1;

  if(tmp){
    return popcount(h) + popcount(d1) + popcount(d2) + popcount(v);
  }
  return 0;
}

void play(long long unsigned int move)
{
	currentPosition ^= mask;
	mask |= move;
  currentPosition |= move;
}
int canPlay(int col)
{
	return (col >= 0 && col < P4_COLONNES && (mask & top_mask(col)) == 0);
}

static long long unsigned int top_mask(int col) 
{
  return (INT64_C(1) << (HEIGHT - 1)) << col*(HEIGHT+1);
}

void playCol(int col)
{
	play((mask + bottom_mask_col(col)) & column_mask(col));
}

void change_current_pos(){
  currentPosition = mask - currentPosition;
}

long long unsigned int alignment(long long unsigned int pos, int n) {
  // horizontal 
  long long unsigned int m = pos & (pos >> (HEIGHT+1));
  if(m = (m & (m >> ((n - 2)*(HEIGHT+1))))) return m;

  // diagonal 1
  m = pos & (pos >> HEIGHT);
  if(m = (m & (m >> ((n - 2)*HEIGHT)))) return m;

  // diagonal 2 
  m = pos & (pos >> (HEIGHT+2));
  if(m = (m & (m >> ((n - 2)*(HEIGHT+2))))) return m;

  // vertical;
  m = pos & (pos >> 1);
  if(m = (m & (m >> (n - 2)))) return m;

  return 0;
}

static long long unsigned int bottom_mask_col(int col) 
{
  return INT64_C(1) << col*(HEIGHT+1);
}

static long long unsigned int column_mask(int col) {
  return ((INT64_C(1) << HEIGHT)-1) << col*(HEIGHT+1);
}

unsigned int entry_index(long long unsigned int key) {
  return (key % 8388593);
}

void put(long long unsigned int key, char val) {
  unsigned int i = entry_index(key); // compute the index position
  entries[i].key = key;              // and overide any existing value.
  entries[i].val = val;       
}

char get(long long unsigned int key) {
  unsigned int i = entry_index(key);  // compute the index position
  if(entries[i].key == key) 
    return entries[i].val;            // and return value if key matches
  else 
    return 0;                   // or 0 if missing entry
}
void reset() { // fill everything with 0, because 0 value means missing data
  memset(&entries[0], 0, 8388593*sizeof(struct Entry));
}

int canWinNext()
{
  return winning_position() & possible();
}

long long unsigned int possible()
{
  return (mask + bottom_mask) & board_mask;
}

int isWinningMove(int col)
{
  return winning_position() & possible() & column_mask(col);
}

long long unsigned int winning_position() {
  return compute_winning_position(currentPosition);
}

long long unsigned int opponent_winning_position() {
  return compute_winning_position(currentPosition /*^ mask*/);
}

static long long unsigned int compute_winning_position(long long unsigned int position) {
  // vertical;
  long long unsigned int r = (position << 1) & (position << 2) & (position << 3);

  //horizontal
  long long unsigned int p = (position << (HEIGHT+1)) & (position << 2*(HEIGHT+1));
  r |= p & (position << 3*(HEIGHT+1));
  r |= p & (position >> (HEIGHT+1));
  p = (position >> (HEIGHT+1)) & (position >> 2*(HEIGHT+1));
  r |= p & (position << (HEIGHT+1));
  r |= p & (position >> 3*(HEIGHT+1));

  //diagonal 1
  p = (position << HEIGHT) & (position << 2*HEIGHT);
  r |= p & (position << 3*HEIGHT);
  r |= p & (position >> HEIGHT);
  p = (position >> HEIGHT) & (position >> 2*HEIGHT);
  r |= p & (position << HEIGHT);
  r |= p & (position >> 3*HEIGHT);

  //diagonal 2
  p = (position << (HEIGHT+2)) & (position << 2*(HEIGHT+2));
  r |= p & (position << 3*(HEIGHT+2));
  r |= p & (position >> (HEIGHT+2));
  p = (position >> (HEIGHT+2)) & (position >> 2*(HEIGHT+2));
  r |= p & (position << (HEIGHT+2));
  r |= p & (position >> 3*(HEIGHT+2));

  return r & (board_mask ^ mask);
}

static long long unsigned int compute_alignment_three_position(long long unsigned int position){
  //vertical
  long long unsigned int r = (position << 1) & (position << 2);

  //horizontal
  r |= (position << (HEIGHT+1)) & (position << 2*(HEIGHT+1));
  r |= (position << (HEIGHT+1)) & (position >> (HEIGHT+1));
  r |= (position >> (HEIGHT+1)) & (position >> 2*(HEIGHT+1));

  //diagonal 1
  r |= (position << HEIGHT) & (position << 2*HEIGHT);
  r |= (position << HEIGHT) & (position >> HEIGHT);
  r |= (position >> HEIGHT) & (position >> 2*HEIGHT);

  //diagonal 2
  r |= (position << (HEIGHT+2)) & (position << 2*(HEIGHT+2));
  r |= (position << (HEIGHT+2)) & (position >> (HEIGHT+2));
  r |= (position >> (HEIGHT+2)) & (position >> 2*(HEIGHT+2));

  return r & (board_mask ^ mask);
}

long long unsigned int possibleNonLoosingMoves() {
  long long unsigned int possible_mask = possible();
  long long unsigned int opponent_win = opponent_winning_position();
  long long unsigned int forced_moves = possible_mask & opponent_win;
  if(forced_moves) {
    if(forced_moves & (forced_moves - 1)) // check if there is more than one forced move
      return 0;                           // the opponnent has two winning moves and you cannot stop him
    else possible_mask = forced_moves;    // enforce to play the single forced move
  }
  return possible_mask & ~(opponent_win >> 1);  // avoid to play below an opponent winning spot
}

static long long unsigned int bottom(int width, int height) {
  return width == 0 ? 0 : bottom(width-1, height) | 1LL << (width-1)*(height+1);
}

int move_score(long long unsigned int move) {
  int count;
  long long unsigned int newPos = ((currentPosition ^ mask) | move);
  return ((count = 50 * popcount(compute_winning_position(newPos))) ? count : popcount(compute_alignment_three_position(newPos)));
}

void add(long long unsigned int move, int score, struct MoveScore* moveScore, int* size)
{
  int pos = (*size)++;
  for(; pos && moveScore[pos-1].score > score; --pos) moveScore[pos] = moveScore[pos-1];
  moveScore[pos].move = move;
  moveScore[pos].score = score;
}

long long unsigned int getNext(struct MoveScore* moveScore, int* size) 
{
  if((*size)) 
    return moveScore[--(*size)].move;
  else 
    return 0;
}

static unsigned int popcount(long long unsigned int m) {
  unsigned int c = 0; 
  for (c = 0; m; c++) m &= m - 1;
  return c;
}

int move_score_end(long long unsigned int move) {
  return popcount(compute_winning_position((currentPosition ^ mask) | move));
}