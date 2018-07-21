#[derive(Debug, Clone)]
pub enum SurfelRule {
    Deteriorate {
        substance_idx: usize,
        factor: f32,
    },
    Transfer {
        source_substance_idx: usize,
        target_substance_idx: usize,
        factor: f32,
    },
}
