/* RCS_INFO = "$RCSfile: read_lp_state.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading LP Files                       */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "iqsutil.h"
#include "read_lp.h"
#include "lp.h"
#include "rawlp.h"
#include "lpdefs.h"
#include "format.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
static int TRACE = 0; 

#define END_LINE(p)  (((*p) == '\\' || (*p) == '\n' || (*p) == '\0') ? 1 : 0)

static const char* all_keyword[] = {
	"MIN", "MINIMUM", "MINIMIZE", 
	"MAX", "MAXIMUM", "MAXIMIZE", 
        "SUBJECT", "ST", "PROBLEM", "PROB",       
	"BOUNDS",  "BOUND", "INTEGER", "END", NULL }; 
static int all_keyword_len[] = {
        3, 7, 8, 
        3, 7, 8, 
	7, 2, 7, 4, 
	6, 5,  7, 3, -1 }; 

int ILLread_lp_state_init (ILLread_lp_state  *state, 
			   qsline_reader *file, const char* fname, 
			   int inter) 
{
    int rval = 0;
    ILL_FAILtrue(file == NULL, "need a file"); 
    state->eof = 0; 
    state->file_name= fname;
    state->interactive= inter; 
    state->file = file; 
    state->line_num = 0; 
    state->p = state->line;  
    state->line[0] = '\0';
    state->realline[0] = '\0';
    state->field[0] = '\0';
    state->fieldOnFirstCol = 0;
    ILLread_lp_state_skip_blanks(state, 1);
CLEANUP:
    ILL_RETURN(rval, "ILLread_lp_state_init"); 
} 

int ILLread_lp_state_next_line (ILLread_lp_state  *state) 
{
    char *slash; 
    if (state->eof) {
	return 1;
    } 
    state->line[0] = '\0';
    if (state->interactive) {
	fprintf(stdout, "> "); fflush(stdout); 
    } 
    while (ILLline_reader_get(state->realline, ILL_namebufsize - 2, state->file) 
	   != (char *) NULL) {
        state->p = state->line; 
        state->line_num++; 
        strcpy(state->line, state->realline); 
        slash = strchr(state->line, '\\'); 
        if (slash != NULL) {
	    *slash = '\0'; 
        } 
        while (ILL_ISBLANK(state->p)) { 
	    state->p++; 
        }
        if (!END_LINE(state->p)) {
	    ILL_IFTRACE("NEWLINE %s %d: %s", 
			state->file_name, state->line_num, state->line); 
	    return 0; 
	}
	if (state->interactive) {
	    fprintf(stdout, "> "); fflush(stdout); 
	} 
    }
    state->eof = 1; 
    state->line_num++; 
    state->field[0] = '\0';
    state->line[0] = '\0';
    strcpy(state->realline, "\n"); 
    state->p = state->line; 
    state->fieldOnFirstCol = 0;
    return 1;
}

int ILLread_lp_state_skip_blanks (ILLread_lp_state  *state, char wrapLines) 
{
    while (1) { 
        while (ILL_ISBLANK(state->p)) { 
	    state->p++; 
        }
        if (END_LINE(state->p)) {
            if (wrapLines) { 
		if (ILLread_lp_state_next_line (state) != 0) {
		    return 1;
		} 
	    } else {
               return 0;  /* done */
            } 
        } else {
            return 0; /* foud non blank */
        } 
    }
}

static int next_field (ILLread_lp_state  *state, int acrossLines)  
{
    (void) ILLread_lp_state_skip_blanks (state, (char) acrossLines); 
    if (state->eof) {
	return 1; 
    } 
    state->fieldOnFirstCol = (state->line == state->p); 
    if (sscanf(state->p, "%s", state->field) != EOF) {
	state->p += strlen(state->field); 
	return 0; 
    }
    return 1; 
}
int ILLread_lp_state_next_field_on_line (ILLread_lp_state  *state)  
{
	return  next_field(state, 0);
} 
int ILLread_lp_state_next_field(ILLread_lp_state  *state)  
{
	return  next_field(state, 1); 
} 

void ILLread_lp_state_prev_field(ILLread_lp_state  *state)  
{
    if (state->p > state->line) {
	state->p--; 
    }  
    while (ILL_ISBLANK(state->p) && (state->p > state->line)) {
        state->p--; 
    } 
    while (!ILL_ISBLANK(state->p) && (state->p > state->line)) {
        state->p--; 
    }
    state->fieldOnFirstCol = (state->line == state->p); 
}

int  ILLread_lp_state_next_var (ILLread_lp_state  *state) 
{
    char *p; 
    int var_len, i; 

    if (ILLread_lp_state_skip_blanks (state, 1)) {
	return 1; 
    } 
    state->fieldOnFirstCol = (state->line == state->p); 
    var_len = 0;
    p = state->p; 
    while (1) { 
	if (ILLis_lp_name_char(*p, var_len)) { 
	    p++; 
            var_len++; 
	} else {
	    break;
	} 
    } 
    if (var_len == 0) { 
        return 1; 
    }
    if (state->fieldOnFirstCol) { 
	/* see whether we founbd a reserved keyword */ 
	for (i= 0; all_keyword[i] != NULL; i++) {
	    if ((var_len == all_keyword_len[i]) && 
		(strncasecmp(all_keyword[i], 
			     state->p, all_keyword_len[i]) == 0)) { 
		return -1;  /* yes we did */
	    } 
	}
    }
    strncpy(state->field, state->p, var_len); 
    state->field[var_len] = '\0'; 
    state->p = p; 
    return 0; 
}

int ILLread_lp_state_bad_keyword (ILLread_lp_state  *state) 
{
    if (!state->fieldOnFirstCol) {
	return  ILLlp_error (state, 
			  "Keyword \"%s\" not at beginning of line.\n", 
			  state->field); 
    }
    return 0; 
}

int ILLtest_lp_state_keyword (ILLread_lp_state  *state, const char* kwd[]) 
{
    int i = 0; 
    if (!state->eof && state->fieldOnFirstCol) { 
	for (i =0; kwd[i] != NULL; i++) {
	    if (strcasecmp(state->field, kwd[i]) == 0) {
		return 0; 
	    }
	}
    }
    return 1; 
} 
int ILLread_lp_state_keyword (ILLread_lp_state  *state, const char* kwd[]) 
{
    if (state->eof || ILLread_lp_state_bad_keyword (state)) {
	return 1; 
    } 
    return  ILLtest_lp_state_keyword (state, kwd); 
} 


int ILLread_lp_state_colon(ILLread_lp_state  *state) 
{
    if ((ILLread_lp_state_skip_blanks (state, 1) == 0) && (*state->p == ':')) {
	state->p++; 
	return 0; 
    } 
    return 1;
}

int  ILLread_lp_state_has_colon (ILLread_lp_state  *state) {
    char *pp; 
    ILLread_lp_state_skip_blanks (state, 0); 
    for (pp = state->p; *pp != '\n'; pp++) {
	if (*pp == ':') {
	    return 1; 
	} 
    } 
    return 0; 
} 

int ILLread_lp_state_next_constraint (ILLread_lp_state  *state) 
{
    int rval; 
    int ln = state->line_num; 

    ILLread_lp_state_skip_blanks (state, 1); 
    if (state->eof) {
        return 1; 
    } 
    if (ln == state->line_num) {
        return ILLlp_error (state, "Constraints must start on a new line.\n"); 
    } 
    if (ILLread_lp_state_next_field (state) == 0) {
       rval = ILLtest_lp_state_keyword (state, all_keyword); 
       ILLread_lp_state_prev_field(state); 
       return !rval;
    } 
    return 0;
} 

/* return 0 if there is a sign */
int ILLread_lp_state_sign (ILLread_lp_state  *state, double *sign)
{
    char found = 0; 
    *sign = 1.0; 
    if (ILLread_lp_state_skip_blanks (state, 1) == 0) {
	if ((*state->p == '+') || (*state->p == '-'))  { 
	    *sign = (*state->p == '+') ? 1.0 : -1.0;
	    state->p++; 
	    found = 1; 
	}
    }
    return 1 - found; 
} 

int ILLtest_lp_state_next_is (ILLread_lp_state  *state, const char *str) 
{
    ILLread_lp_state_skip_blanks (state, 0); 
    if (strncasecmp(state->p, str, strlen(str)) == 0) {
        state->p += strlen(str); 
	return 1; 
    }
    return  0; 
} 

int ILLread_lp_state_value (ILLread_lp_state  *state, double *coef) 
{
    int len = 0; 

    if (ILLread_lp_state_skip_blanks (state, 1) != 0) {
        ILL_RESULT(1, "ILLread_lp_state_value"); 
    } else {
        state->fieldOnFirstCol = (state->line == state->p); 
        len = ILLget_value(state->p, coef); 
        if (len > 0) {
            state->p += len; 
            ILL_RESULT(0, "ILLread_lp_state_value"); 
        }
        ILL_RESULT(1, "ILLread_lp_state_value"); 
    }
}

int  ILLread_lp_state_possible_coef (ILLread_lp_state  *state, double *coef, 
                                        double defValue) 
{
    *coef = defValue; 
    return  ILLread_lp_state_value (state, coef);  
} 


int  ILLread_lp_state_possible_bound_value (ILLread_lp_state  *state) 
{
    double sign; 
    int len = 0; 
    char *p = NULL; 
    (void) ILLread_lp_state_sign (state, &sign); 

    if (!strncasecmp (state->p, "INFINITY", 8)) { 
	len = 8; 
    } else { 
	if (!strncasecmp (state->p, "INF", 3)) {
	    len = 3;
	}
    }
    if (len > 0) {
        state->p += len; 
        p = state->p;
        ILLread_lp_state_skip_blanks (state, 0); 
        if (!END_LINE(p) && p == state->p) { 
            /* found no blanks so this INF/INFINITY is the prefix 
             * of something else */
            state->p -= len; 
            return 0; /* no coef found */ 
        } else {
	    state->bound_val = sign * ILL_MAXDOUBLE;
	    return 1; 
        }
    } 
    if (ILLread_lp_state_value (state, &state->bound_val) == 0) {
	state->bound_val *= sign; 
	return 1; 
    }
    return 0; /* no coef found */ 
}

int ILLtest_lp_state_sense (ILLread_lp_state  *state, char all) {
    char c; 

    state->sense_val = ' '; 
    if (ILLread_lp_state_skip_blanks (state, 1) == 0) {
	c = *state->p; 
	if (!all) {   /* look for '=' and '<=' */
	    if (c == '=') { 
		state->p++;
		state->sense_val = 'E'; 
	    } else { 
                if ((c == '<') && (*(state->p+1) == '=')) {
		    state->p += 2; 
		    state->sense_val = 'L'; 
		} 
	    }
        } else { 
	    c = *state->p; 
	    if ((c == '<') || (c == '>')) {
		state->sense_val = (c == '<') ? 'L' : 'G';
		state->p++; 
		c = *state->p; 
		if (*state->p == '=') {
		    state->p++; 
		} 
	    } else { 
                if (c == '=') { 
		    state->p++; 
		    c = *state->p; 
		    if ((c == '<') || (c == '>')) {
			state->sense_val = (c == '<') ? 'L' : 'G';
			state->p++; 
		    } else { 
			state->sense_val = 'E';
		    }
		} 
	    }
	}
    }
    return (state->sense_val != ' '); 
}

void ILLtest_lp_state_bound_sense (ILLread_lp_state  *state)
{
    (void) ILLtest_lp_state_sense (state, 0); 
}

int ILLread_lp_state_sense (ILLread_lp_state  *state)
{
    if (! ILLtest_lp_state_sense (state, 1)) {
	if (END_LINE(state->p)) { 
	    return ILLlp_error (state, "Missing row sense at end of line.\n");
	} else { 
            if (*state->p != '\0') { 
		return ILLlp_error (state, "\"%c\" is not a row sense.\n", 
				    *state->p);
            } else {
		return ILLlp_error (state, "Missing row sense at end of line.\n");
	    }
	}
    } 
    return 0; 
}

/* ------------------------------------------------------------------------- */
/*  error printing 
 */

static void ILLread_lp_state_print_at(ILLread_lp_state  *state) {
    char *p; 
    if (state->eof) {
	fprintf (stderr, "end of file"); 
    } else { 
        if (*state->p == '\n') {
	    fprintf (stderr, "end of line"); 
	}  else {
	    p = state->p; 
	    while (ILL_ISBLANK(p)) {  
		p++; 
	    }
	    fprintf(stderr, "%c", '"'); 
	    for (; !ILL_ISBLANK(p) && !END_LINE(p); p++) {  
		fprintf(stderr, "%c", *p); 
	    }
	    fprintf(stderr, "\""); 
	} 
    }
}

static void lp_err(ILLread_lp_state *state, int isError, 
		   const char* format, va_list args)
{ 
    int rval = 0;
    int errtype, slen, at;
    qsformat_error error; 
    char error_desc[256]; 

    ILL_FAILfalse(state != NULL, "state != NULL"); 
    ILL_FAILfalse(state->file != NULL, "state->file != NULL"); 
    ILL_FAILfalse(format != NULL, "format != NULL"); 
    ILL_FAILfalse(format[0] != '\0', "format[0] != '0'"); 

    ILLread_lp_state_skip_blanks (state, 0); 
    at = state->p - state->line; 
    vsprintf(error_desc, format, args); 
    slen = strlen(error_desc); 
    if ((slen > 0) && error_desc[slen-1] != '\n') {
	error_desc[slen] = '\n'; 
	error_desc[slen+1] = '\0'; 
    }

    if (state->file->error_collector != NULL) {
	errtype = (isError) ? QS_LP_FORMAT_ERROR : QS_LP_FORMAT_WARN;
	ILLformat_error_create(&error, errtype, error_desc, 
			       state->line_num, state->realline,  at);
        ILLformat_error(state->file->error_collector, &error); 
        ILLformat_error_delete(&error); 
    } else { 
	if (!state->interactive){
	    fprintf(stderr, "%s %d: %s\t", state->file_name, state->line_num, 
		    state->realline); 
	    fprintf (stderr, "%s at ", 
		     (isError) ? "LP Error" : "LP Warning"); 
	    ILLread_lp_state_print_at(state); 
	    fprintf (stderr, ": "); 
	} else {
	    fprintf (stderr, "%s : ",  (isError) ? "LP Error" : "LP Warning");
	} 
        fprintf(stderr, error_desc); 
        fflush(stderr); 
    }
CLEANUP: ;
}

int ILLlp_error (ILLread_lp_state  *state, const char* format, ...) 
{
    va_list args;
    va_start(args, format);
    lp_err(state, TRUE,format, args); 
    return 1; 
}

void ILLlp_warn (ILLread_lp_state  *state, const char* format, ...) 
{
    va_list args;
    va_start(args, format);
    if (format != NULL) { 
	lp_err(state, FALSE, format, args); 
    }
}

/* shared with read_mps_state.c */
int ILLget_value (char *line, double *coef) 
{
	char field[ILL_namebufsize];
	int rval = 0, i;
	char c, lastC, *p;
	int allowDot, allowExp, allowSign;

	p = line;
	c = *p;
	i = 0;
	lastC = ' '; 
	allowSign = 1;
	allowExp = 0;
	allowDot = 1;
	while ((('0' <= c) && (c <= '9')) ||
		(allowDot && (c == '.')) || 
		((allowExp == 1) && ((c == 'e') || (c == 'E'))) ||
		((allowSign || lastC == 'e' || lastC == 'E') && 
		((c == '+') || (c == '-'))) ) {
			if (c == '.') allowDot = 0;
			allowSign = 0;

			if ((allowExp == 0) && (c >= '0') && (c <= '9')) {
				allowExp = 1; 
			}
			if ((c == 'e') || (c == 'E')) {
				allowExp++; 
				allowDot = 0;
			} 
			p++;
			lastC = c; 
			c = *p;
			i++;
		}
		if ((lastC == '+') || (lastC == '-')) {
			p--; i--; 
			if (p > line) 
				lastC = *(p-1);
			else 
				lastC = ' '; 
		}
		if ((lastC == 'e') || (lastC == 'E')) {
			p--;
			i--;
		}
		if (i > 0) {
			strncpy(field, line, i);
			field[i] = '\0';
			rval = !sscanf(field, "%lf%n", coef, &i); 
			if (rval != 0) {
				ILL_RESULT(0, "ILLget_value"); 
			}
		}
		ILL_RESULT(i, "ILLget_value"); 
}

int ILLcheck_subject_to(ILLread_lp_state *state) 
{
    int rval; 
    char *p; 
    if ((rval = ILLread_lp_state_next_field(state)) == 0) {
	if (strcasecmp (state->field, "ST") == 0) { 
	    rval = ILLread_lp_state_bad_keyword(state); 
	}  else { 
	    if (strcasecmp (state->field, "SUBJECT") == 0) { 
		p = state->p; 
		while (ILL_ISBLANK(p)) { 
		    p++; 
		}
		if (!strncasecmp (p, "TO", 2)) { 
		    rval =  ILLread_lp_state_bad_keyword(state); 
		    if (rval == 0) {
			state->p = p + 2; 
		    } 
		}
	    } else {
		rval = 1; 
	    } 
	}
	if (rval != 0) {
	    ILLread_lp_state_prev_field(state); 
	}  else { 
            ILLread_lp_state_skip_blanks(state, 1); 
        }
    }
    ILL_RESULT(rval, "check_subject_to"); 
}

