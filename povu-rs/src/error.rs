//! Error types for the Povu Rust API

use std::ffi::NulError;
use std::fmt;

/// Result type alias for Povu operations
pub type Result<T> = std::result::Result<T, Error>;

/// Error type for Povu operations
#[derive(Debug)]
pub enum Error {
    /// Error from the underlying C++ library
    Povu { code: i32, message: String },
    /// Null pointer or invalid data
    NullPointer,
    /// String conversion error
    StringConversion(String),
    /// I/O error
    Io(std::io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Povu { code, message } => {
                write!(f, "Povu error (code {}): {}", code, message)
            }
            Error::NullPointer => write!(f, "Null pointer encountered"),
            Error::StringConversion(msg) => write!(f, "String conversion error: {}", msg),
            Error::Io(e) => write!(f, "I/O error: {}", e),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<NulError> for Error {
    fn from(e: NulError) -> Self {
        Error::StringConversion(format!("Null byte in string: {}", e))
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

impl Error {
    /// Create an Error from a Povu FFI error
    pub(crate) fn from_ffi_error(error: crate::ffi::PovuError) -> Self {
        if error.message.is_null() {
            Error::Povu {
                code: error.code,
                message: "Unknown error".to_string(),
            }
        } else {
            let message = unsafe {
                std::ffi::CStr::from_ptr(error.message)
                    .to_string_lossy()
                    .into_owned()
            };
            // Free the error message
            unsafe {
                crate::ffi::povu_string_free(error.message);
            }
            Error::Povu {
                code: error.code,
                message,
            }
        }
    }
}
