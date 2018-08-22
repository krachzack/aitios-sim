use self::Transport::*;
use std::default::Default;
use transport::{transport, Absorb, AbsorbThenDeposit, Deposit, DepositAll, Transport as Inner};

/// Specifies when and in which direction substance is transported.
pub enum Transport {
    /// Absorb on bounce, deposit on settle.
    Classic(Inner<Absorb, Deposit>),
    /// On both bounce and settle, absorb first and then deposit.
    Consistent(Inner<AbsorbThenDeposit, AbsorbThenDeposit>),
    /// Disposes of all substance on settle, otherwise consistent
    Conserving(Inner<AbsorbThenDeposit, DepositAll>),
}

impl Transport {
    pub fn classic() -> Self {
        Classic(transport())
    }

    pub fn consistent() -> Self {
        Consistent(transport())
    }

    pub fn conserving() -> Self {
        Conserving(transport())
    }
}

impl Default for Transport {
    fn default() -> Self {
        Self::classic()
    }
}
