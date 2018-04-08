use surfel_rule::SurfelRule;

#[derive(Debug, Clone)]
pub struct SurfelData {
    /// Index of the entity in the scene that the triangle that this surfel was generated from belongs to.
    /// Important to assign surfels to textures.
    pub entity_idx: usize,
    /// Deterioration rate of the probability of a gammaton moving further in a straight line.
    pub delta_straight: f32,
    /// Deterioration rate of the probability of a gammaton moving in a piecewise approximated parabolic path.
    pub delta_parabolic: f32,
    /// Deterioration rate of the probability of a gammaton flowing in a tangent direction.
    pub delta_flow: f32,
    /// Holds the amount of substances as numbers in the interval 0 ≤ a ≤ 1.
    pub substances: Vec<f32>,
    /// Weights for the transport of substances from a settled ton to a surfel.
    pub deposition_rates: Vec<f32>,
    pub rules: Vec<SurfelRule>
}
