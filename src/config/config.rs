use config::Transport;

/// Encapsulates parameters that influence substance transport and tracing.
#[derive(Default)]
pub struct Config {
    pub transport: Transport,
}
