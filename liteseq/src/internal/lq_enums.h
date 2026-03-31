#ifndef LQ_ENUMS_H
#define LQ_ENUMS_H

#define DEFINE_ENUM(enum_name, ENUM_DEF)                                       \
	enum enum_name { ENUM_DEF(ENUM_ITEM) enum_name##_INVALID };            \
	/* Function prototypes */                                              \
	const char *to_string_##enum_name(enum enum_name value);               \
	enum enum_name from_string_##enum_name(const char *s);

#define DEFINE_ENUM_AND_STRING(enum_name, ENUM_DEF)                            \
	/* Function definitions */                                             \
	const char *to_string_##enum_name(enum enum_name value)                \
	{                                                                      \
		switch (value) {                                               \
			ENUM_DEF(ENUM_CASE) /* Maps enum -> string */          \
		default:                                                       \
			return "UNKNOWN_" #enum_name;                          \
		}                                                              \
	}                                                                      \
									       \
	enum enum_name from_string_##enum_name(const char *s)                  \
	{                                                                      \
		if (s) {                                                       \
			ENUM_DEF(ENUM_COMPARE) /* Maps string -> enum */       \
		}                                                              \
		return enum_name##_INVALID;                                    \
	}

/*
 * sym     is the enum constant (e.g., GFA_1_0)
 * literal is the paired string (e.g., "VN:Z:1.0")
 * s       is the runtime query string passed into from_string_*
 */
#define ENUM_ITEM(sym, literal) sym, // Enum generation
#define ENUM_CASE(sym, literal)                                                \
	case sym:                                                              \
		return literal; // For to_string
#define ENUM_COMPARE(sym, literal)                                             \
	if (strcmp(s, literal) == 0)                                           \
		return sym; // For from_string

#endif /* LQ_ENUMS_H */
