#pragma once

#ifndef _PI 
#define	_PI (double)(4.0*atan((double)1))
#endif

#define _N 0x1000


/********** User IDs for callbacks ********/
#define LIGHT0_ENABLED_ID    200
#define LIGHT1_ENABLED_ID    201
#define LIGHT0_INTENSITY_ID  250
#define LIGHT1_INTENSITY_ID  260
#define ENABLE_ID            300
#define DISABLE_ID           301
#define HIDE_ID              303

#define CK_EVENLY_SPACED          1000
#define CK_WINDOW_INVERSE         1001
#define BT_NOISE_DEFAULT          1002
#define BT_SPHERE_ADJUST          1003
#define ET_MERIDIAN_MAIN_MATERIAL 1004
#define ET_MERIDIAN_MINOR_MATERIAL  1005
#define ET_PARALLEL_MAIN_MATERIAL 1006
#define ET_PARALLEL_MINOR_MATERIAL  1007
#define ET_DIAGONALE_UR_MATERIAL  1008
#define ET_DIAGONALE_DR_MATERIAL  1009
#define ET_SPHERE_MP_MATERIAL     1010
#define ET_SPHERE_Mp_MATERIAL     1011
#define ET_SPHERE_mP_MATERIAL     1012
#define ET_SPHERE_mp_MATERIAL     1013
#define CK_FACE                   1014
#define CK_SMOOTH                 1015
#define CK_WIRE                   1016
#define ET_FACE_CENTER_MATERIAL   1017
#define ET_DOME_NAME              1018
#define BT_DOME_LOAD              1019
#define BT_DOME_SAVE              1020
#define BT_DOME_INC               1021
#define ET_DOME_MATERIAL          1022
#define RG_DISPLAY                1023
#define CK_DISPLAY_MINOR          1024
#define CK_LIGHT_1                1025
#define CK_LIGHT_2                1026
#define BT_DEFAULT_VIEW           1027
#define BT_TOP_VIEW               1028
#define BT_SIDE_VIEW              1029
#define SP_WINDOW_PERSIDE         1030
#define LB_BASE                   1031
#define SP_BASE_SIDE			  1032
#define SP_BASE_PARAMETER		  1033
#define SP_BASE_SMOOTHING         1034
#define LB_SLOPE				  1035
#define SP_SLOPE_TWIST			  1036
#define SP_DOME_SCALE_Z			  1037
#define SP_MERIDIAN_MAIN		  1038
#define	SP_MERIDIAN_MINOR   	  1039
#define	SP_PARALLEL_MAIN		  1040
#define	SP_PARALLEL_MINOR   	  1041
#define LB_WINDOW				  1042
#define SP_WINDOW_HEIGHT		  1043
#define SP_WINDOW_SHIFT			  1044
#define CK_DISPLAY_SMOOTH         1045
#define SP_NOISE_OCTAVE           1046
#define SP_NOISE_DENSITY          1047
#define SP_NOISE_SEED             1048
#define SP_NOISE_FREQUENCY_U      1049
#define SP_NOISE_FREQUENCY_V      1050
#define SP_NOISE_AMPLITUDE_R      1051
#define SP_NOISE_AMPLITUDE_Z      1052
#define SP_NOISE_PROGRESSIVE_R    1053
#define SP_NOISE_PROGRESSIVE_Z    1054
#define RG_SLOPE_SPACED           1056
#define SP_WINDOW_DELAY           1057
#define SP_WINDOW_ATTACK          1058
#define SP_DIAGONAL_UR			  1059
#define SP_DIAGONAL_DR			  1060
#define SP_WINDOW_DAMPING         1061
#define SP_BASE_EXPAND            1062
#define LB_MATERIAL               1063
#define BT_QUIT_MESSAGE           1064
#define SP_SLOPE_PARAMETER        1065
#define BT_DOME_DEFAULT           1066
#define BT_SLOPE_DEFAULT          1067
#define BT_BASE_DEFAULT           1068
#define BT_WINDOW_DEFAULT         1069
#define SP_BASE_DAMPING           1070
#define SP_SLOPE_PARAMETER_XY	  1071
#define SP_WINDOW_POWER  		  1072
#define BT_DOME_STL               1073
#define BT_DOME_OBJ               1074
#define SP_DOME_END				  1075
#define SP_DOME_SOCLE             1076
#define LB_DOME_FILE			  1077
#define CK_LISSAGE				  1078
#define SP_LISSAGE_POIDS_SELF	  1079
#define SP_LISSAGE_POIDS_HORIZONTAL  1080
#define SP_LISSAGE_POIDS_DIAGONALES		  1081
#define SP_LISSAGE_NOMBRE		  1082
#define SP_LISSAGE_TAILLE		  1083
#define BT_LISSAGE_GAUSSIEN	      1084
#define BT_LISSAGE_LAPLACIEN	  1085
#define BT_LISSAGE_DEFAULT  	  1086
#define SP_LISSAGE_DISTANCE_TYPE  1087
#define SP_LISSAGE_PUISSANCE	  1088
#define CK_EVOLUTE				  1089
#define SP_EVOLUTE_LAMBDA		  1090
#define SP_EVOLUTE_NOMBRE         1091
#define CK_EVOLUTE_MERIDIEN       1092
#define CK_EVOLUTE_PARALLEL       1093
#define CK_NORMALE				  1094
#define SP_NORMALE_LAMBDA		  1095
#define SP_NORMALE_NOMBRE         1096
#define SP_SPIRALE				  1097
#define LB_SPIRALE_SHRINK		  1098
#define SP_SPIRALE_RADIUS		  1099
#define ET_SPIRALE_MATERIAL		  1100
#define SP_SPIRALE_TOUR			  1101
#define SP_SPIRALE_SPEED		  1102
#define SP_SPIRALE_DELAY		  1103
#define SP_SPIRALE_DECALAGE		  1104
#define LB_SPIRALE_CHIRALE		  1105
#define SP_FRISE				  1106
#define SP_FRISE_HEIGHT			  1107
#define SP_FRISE_PERSIDE		  1108
#define SP_FRISE_START			  1109
#define SP_FRISE_SHIFT			  1110
#define SP_FRISE_RADIUS			  1111
#define LB_FRISE_SHRINK			  1112
#define LB_FRISE_FUNCTION		  1113
#define ET_FRISE_MATERIAL		  1114

#define SP_SNAKE				  1115
#define SP_SNAKE_TOUR			  1116		
#define LB_SNAKE_EVOLVE		      1117
#define SP_SNAKE_SHIFT			  1118
#define SP_SNAKE_RADIUS			  1119
#define LB_SNAKE_SHRINK			  1120
#define ET_SNAKE_MATERIAL		  1121
#define LB_SNAKE_CHIRALE		  1122
#define SP_SNAKE_STEP			  1123
#define SP_SNAKE_OFFSET			  1124
#define LB_SNAKE_TOWARD			  1125

#define ET_FACE_EDGE_MATERIAL     1126
#define BT_OK_MESSAGE             1127

#define SP_MERIDIAN_MAIN_RADIUS   1128
#define SP_MERIDIAN_MINOR_RADIUS  1129
#define SP_PARALLEL_MAIN_RADIUS   1130
#define SP_PARALLEL_MINOR_RADIUS  1131
#define SP_SPHERE_MP_RADIUS       1132
#define SP_SPHERE_Mp_RADIUS       1133
#define SP_SPHERE_mP_RADIUS       1134
#define SP_SPHERE_mp_RADIUS       1135
#define SP_WIRE_SHRINKING         1136
#define SP_WIRE_SHRINKING_SPEED   1137
#define LB_WIRE_SHRINKING_TYPE    1138

#define SP_MERIDIAN_EDGE_SIZE     1139
#define SP_PARALLEL_EDGE_SIZE     1140

#define SP_SPIRALE_STEP			  1141
#define LB_SPIRALE_TOWARD		  1142

#define LB_FRISE_TOWARD			  1143
#define SP_FRISE_STEP			  1144
#define SP_FRISE_POWER			  1145

#define SP_DERIVATIVE			  1146
#define CK_DERIVATIVE_NORME		  1147

#define BT_DECOR_QUIT			  1148
#define BT_POST_PROCESSING_QUIT   1149
#define BT_DECOR				  1150
#define BT_POST_PROCESSING		  1151

#define SP_FRISE_OFFSET			  1152

#define SP_BUMP_DAMPING			  1153
#define SP_BUMP_POWER			  1154
#define CK_BUMP_NORMALIZE		  1155
#define LB_BUMP_NORMAL			  1156
#define SP_BUMP_LAMBDA			  1157
#define CK_BUMP					  1158
#define SP_BUMP_DISTANCE		  1159
#define LB_BUMP_DISTANCE_TYPE	  1160
#define CK_BUMP_DISTANCE_INVERSE  1161

#define SP_SPIKE_MP				  1162
#define SP_SPIKE_Mp				  1163
#define SP_SPIKE_mP				  1164
#define SP_SPIKE_mp				  1165
#define BT_SPIKE_DEFAULT		  1166

#define BT_PRE_PROCESSING_QUIT	  1167
#define BT_PRE_PROCESSING		  1168

#define ET_LID_MATERIAL			  1169
#define CK_LID					  1170
#define SP_DOME_TOP				  1171
#define CK_SLOPE_INVERSE		  1172

#define BT_BASE                   1173
#define BT_SLOPE				  1174
#define BT_WINDOW				  1175

#define BT_BASE_QUIT              1176
#define BT_SLOPE_QUIT			  1177
#define BT_WINDOW_QUIT			  1178

#define SP_BASE_ROUND			  1179
#define SP_BASE_BARYCENTRE		  1180

#define LB_BASE_2				  1181

#define SP_CRENEL_PERSIDE		  1182
#define SP_CRENEL_FRONT			  1183
#define SP_CRENEL_THICK			  1184
#define SP_CRENEL_WIDE			  1185
#define SP_CRENEL_HEIGHT		  1186
#define SP_CRENEL_POWER			  1187

#define SP_WINDOW_AWNING		  1188
#define SP_WINDOW_PARAMETER	      1189
#define SP_SLOPE_DELAY			  1190
#define SP_BASE_MORPHING		  1191
#define SP_BASE_SECONDARY_ROTATION 1192
#define SP_BASE_SECONDARY_EXPANSION 1193
#define CK_BASE_MORPHING_DIVISION  1194
#define CK_WINDOW_BLIND           1195
#define SP_FRISE_GAP			  1196
#define SP_FRISE_PARAMETER		  1197
#define LB_TEST_SLOPE_XY          1198
#define LB_TEST_SLOPE_Z           1199
#define SP_BASE_SPIRAL			  1200
#define LB_SNAKE_SLOPE			  1201
#define SP_SNAKE_SLOPE_PARAMETER  1202
#define LB_SPIRALE_SLOPE		  1203
#define SP_SPIRALE_SLOPE_PARAMETER  1204
#define LB_LISSAGE_FILTRE		  1205
#define SP_LISSAGE_SIGMA_1		  1206
#define SP_LISSAGE_SIGMA_2		  1207
#define SP_LISSAGE_POIDS_VERTICAL 1208
#define LB_LISSAGE_DISTANCE_TYPE  1209
#define LB_LISSAGE_DIRECTION      1210
#define LB_DISPLAY				  1211
#define LB_TOP_TYPE				  1212
#define SP_TOP_PARAMETER		  1213
#define BT_TOP_QUIT				  1214
#define BT_TOP					  1215
#define BT_ADJUST_SCALE			  1216
#define SP_DOME_INFOLDING		  1217

#define LB_BUMP_TYPE			  1218

#define SP_POST_MERIDIAN_EVERY	  1219
#define SP_POST_MERIDIAN_DELAY	  1220
#define SP_POST_PARALLEL_EVOLVE   1221
#define SP_POST_PARALLEL_FIRST	  1222
#define SP_POST_PARALLEL_EVERY	  1223

#define SP_CRENEL_DELAY			  1224
#define SP_CRENEL_ANGLE			  1225
#define RG_BASE_SPACED			  1226

#define SP_SPIKE_ANGLE			  1227

#define CK_FRISE_INVERSE		  1228
#define CK_FRISE_REVERSE		  1229

#define CK_MARQUEES				  1230
#define LB_MARQUEE_FUNCTION		  1231
#define SP_MARQUEE_PARAMETER	  1232
#define SP_MARQUEE_LAMBDA		  1233
#define SP_MARQUEE_DAMPING		  1234
#define SP_MARQUEE_ALPHA		  1235
#define SP_MARQUEE_BETA			  1236
#define CK_MARQUEE_NORMALIZE	  1237
#define SP_MARQUEE_TOP			  1238
#define CK_POST_MERIDIAN_PER_SIDE 1239
#define LB_DIFF_GEOM			  1240

#define LB_ENHANCE_TYPE			  1241
#define CK_ENHANCE_MERIDIAN		  1242
#define CK_ENHANCE_PARALLEL		  1243
#define LB_ENHANCE_MERGE		  1244
#define LB_ENHANCE_CALCUL		  1245
#define SP_ENHANCE_LAMBDA		  1246
#define SP_ENHANCE_DAMPING		  1247

#define SP_WINDOW_LOWER			  1248

#define SP_COMB					  1249
#define SP_COMB_RADIUS			  1250
#define LB_COMB_SHRINK			  1251
#define ET_COMB_MATERIAL		  1252
#define SP_COMB_SHIFT			  1253
#define SP_COMB_STEP			  1254
#define LB_COMB_TOWARD			  1255

#define SP_FRISE_SPACING		  1257
#define BT_FRISE_ADJUST			  1258

#define LB_COMB_TYPE			  1259
#define SP_COMB_THREAD			  1260
#define SP_COMB_THREAD_GAP		  1261
#define SP_COMB_THREAD_GAP_EVOLVE 1262


#define CK_WIRE_FRAME			  1500
#define CK_FLAT_FACE			  1501
#define CK_SMOOTH_FACE			  1502

#define RB_WIRE                  0
#define RB_FLAT                  1
#define RB_SMOOTH                2

#define RB_ALONG_PARAMETER       0
#define RB_ALONG_LENGTH          1
#define RB_ALONG_HEIGHT          2
#define RB_ALONG_ANGLE			 2

#define S_NONE "None"
#define MAX_COULEUR  40

#define SAMPLE_SIZE 1024
#define PerlinNoiseLerp(t, a, b) ( a + t * (b - a) )
#define PerlinNoiseS_curve(t) ( t * t * (3.0f - 2.0f * t) )
#define PerlinNoiseSetup(i,b0,b1,r0,r1)\
	t = vec[i] + _N;\
	b0 = ((int)t) & (SAMPLE_SIZE-1);\
	b1 = (b0+1) & (SAMPLE_SIZE-1);\
	r0 = t - (int)t;\
	r1 = r0 - 1.0f;
