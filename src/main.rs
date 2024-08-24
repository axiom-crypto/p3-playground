use p3_baby_bear::{BabyBear, DiffusionMatrixBabyBear};
use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractField, Field};
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_matrix::dense::RowMajorMatrix;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_playground::fib_air::air::FibonacciAir;
use p3_playground::fib_air::trace::generate_trace_rows;
use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::{prove, verify, StarkConfig};
use rand::thread_rng;

// STARK setup copied from https://github.com/Plonky3/Plonky3/blob/55832146c86e8e4d246bb9843da17f2159d212a5/uni-stark/tests/fib_air.rs#L104
type Val = BabyBear;
type Perm = Poseidon2<Val, Poseidon2ExternalMatrixGeneral, DiffusionMatrixBabyBear, 16, 7>;
type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
type ValMmcs =
    FieldMerkleTreeMmcs<<Val as Field>::Packing, <Val as Field>::Packing, MyHash, MyCompress, 8>;
type Challenge = BinomialExtensionField<Val, 4>;
type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
type Dft = Radix2DitParallel;
type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;

/// Trace generation for the n-th Fibonacci number, where `x` is the n-th Fibonacci number
fn generate_trace_and_pvs(n: usize, x: usize) -> (RowMajorMatrix<Val>, Vec<Val>) {
    let a = 0;
    let b = 1;
    let pvs = [a, b, x].map(BabyBear::from_canonical_usize).to_vec();
    let trace = generate_trace_rows::<Val>(a as u32, b as u32, n);

    (trace, pvs)
}

fn main() {
    // === SETUP ===
    let perm = Perm::new_from_rng_128(
        Poseidon2ExternalMatrixGeneral,
        DiffusionMatrixBabyBear::default(),
        &mut thread_rng(),
    );
    let hash = MyHash::new(perm.clone());
    let compress = MyCompress::new(perm.clone());
    let val_mmcs = ValMmcs::new(hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let dft = Dft {};
    let fri_config = FriConfig {
        log_blowup: 2,
        num_queries: 28,
        proof_of_work_bits: 8,
        mmcs: challenge_mmcs,
    };
    let pcs = Pcs::new(dft, val_mmcs, fri_config);
    let config = MyConfig::new(pcs);
    let mut challenger = Challenger::new(perm.clone());

    // === TRACE GENERATION ===
    let (trace, pvs) = generate_trace_and_pvs(1 << 3, 21);

    // === PROVE ===
    let proof = prove(&config, &FibonacciAir {}, &mut challenger, trace, &pvs);

    // === VERIFY ===
    let mut challenger = Challenger::new(perm);
    verify(&config, &FibonacciAir {}, &mut challenger, &proof, &pvs).expect("verification failed");
}
