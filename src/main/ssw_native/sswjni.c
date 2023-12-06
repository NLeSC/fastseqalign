
#include <jni.h>
#include <stdio.h>
#include "ssw.h"


jobject s_align_to_ssw_Alignment(JNIEnv* env, s_align* align) {
	if (align == NULL) return NULL;
	jclass clazz = (*env)->FindClass(env, "nl/escience/alignment/Alignment");
	if (clazz == NULL) fprintf(stderr,"class cannot be found\n");
	jmethodID constructor = (*env)->GetMethodID(env, clazz, "<init>", "(SIIII)V");
	jobject result = (*env)->NewObject(env, clazz, constructor,
		align->score,
		align->ref_begin,
		align->ref_end,
		align->read_begin,
		align->read_end);

	return result;
}
JNIEXPORT jobject JNICALL Java_nl_escience_alignment_LocalSequenceAlignment_alignNative(JNIEnv* env, jclass cls,
		jbyteArray read, jbyteArray matrix, jint matrixSize,
		jbyteArray ref,
		jint gapOpen,
		jint gapExtend) {
	jbyte* readPtr = (*env)->GetByteArrayElements(env, read, NULL);
	jsize readLen = (*env)->GetArrayLength(env, read);
	jbyte* matrixPtr = (*env)->GetByteArrayElements(env, matrix, NULL);
	/*jsize matrixLen = (*env)->GetArrayLength(env, matrix);*/
	jbyte* refPtr = (*env)->GetByteArrayElements(env, ref, NULL);
	jsize refLen = (*env)->GetArrayLength(env, ref);
	s_profile* profile = ssw_init(readPtr, readLen, matrixPtr, matrixSize);
	s_align* align = ssw_align(profile, refPtr, refLen, gapOpen, gapExtend);
	jobject jalignment = s_align_to_ssw_Alignment(env, align);
	align_destroy(align);
	init_destroy(profile);
	(*env)->ReleaseByteArrayElements(env, read, readPtr, JNI_ABORT);
	(*env)->ReleaseByteArrayElements(env, matrix, matrixPtr, JNI_ABORT);
	(*env)->ReleaseByteArrayElements(env, ref, refPtr, JNI_ABORT);
	return jalignment;
}