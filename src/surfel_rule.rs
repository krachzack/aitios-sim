#[derive(Debug, Clone)]
pub enum SurfelRule {
    // Add a multiple of the current concentration of a substance after every iteration.
    Deteriorate {
        substance_idx: usize,
        factor: f32,
    },
    /// Translate one substance 1:1 to another substance
    Transfer {
        source_substance_idx: usize,
        target_substance_idx: usize,
        factor: f32,
    },
    /// Add a constant amount of substance after every iteration.
    Deposit {
        substance_idx: usize,
        amount: f32
    },
}
