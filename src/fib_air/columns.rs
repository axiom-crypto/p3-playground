pub const NUM_FIBONACCI_COLS: usize = 2;

#[repr(C)]
pub struct FibonacciCols<F> {
    pub left: F,
    pub right: F,
}

impl<F> FibonacciCols<F> {
    pub const fn new(left: F, right: F) -> FibonacciCols<F> {
        FibonacciCols { left, right }
    }
}
