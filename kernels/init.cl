typedef struct {
    ulong hash;
    short idx;
} hash_entry;

__kernel void init( __global hash_entry * g_htable, __global hash_entry * g_htable_even, __global hash_entry * g_htable_odd, const ulong size ) {

	const uint gid = get_global_id(0);

	if(gid<size){
		g_htable[gid].idx = -1;
		g_htable_even[gid].idx = -1;
		g_htable_odd[gid].idx = -1;
	}
}


