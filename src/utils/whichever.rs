/// Helper macro for [`define_whichever`], responsible for implementing `Read`.
macro_rules! impl_read {
    ($struct_name:ident, $match_macro_name:ident) => {
        impl ::std::io::Read for $struct_name {
            #[inline]
            fn read(&mut self, buf: &mut [u8]) -> ::std::io::Result<usize> {
                $match_macro_name!(self, inner => inner.read(buf))
            }

            #[inline]
            fn read_vectored(&mut self, bufs: &mut [::std::io::IoSliceMut<'_>]) -> ::std::io::Result<usize> {
                $match_macro_name!(self, inner => inner.read_vectored(bufs))
            }

            #[inline]
            fn read_to_end(&mut self, buf: &mut ::std::vec::Vec<u8>) -> ::std::io::Result<usize> {
                $match_macro_name!(self, inner => inner.read_to_end(buf))
            }

            #[inline]
            fn read_to_string(&mut self, buf: &mut ::std::string::String) -> ::std::io::Result<usize> {
                $match_macro_name!(self, inner => inner.read_to_string(buf))
            }

            #[inline]
            fn read_exact(&mut self, buf: &mut [u8]) -> ::std::io::Result<()> {
                $match_macro_name!(self, inner => inner.read_exact(buf))
            }
        }
    }
}

/// Helper macro for [`define_whichever`], responsible for implementing `Write`.
macro_rules! impl_write {
    ($struct_name:ident, $match_macro_name:ident) => {
        impl ::std::io::Write for $struct_name {
            #[inline]
            fn write(&mut self, buf: &[u8]) -> ::std::io::Result<usize> {
                $match_macro_name!(self, inner => inner.write(buf))
            }

            #[inline]
            fn flush(&mut self) -> ::std::io::Result<()> {
                $match_macro_name!(self, inner => inner.flush())
            }

            fn write_vectored(&mut self, bufs: &[::std::io::IoSlice<'_>]) -> ::std::io::Result<usize> {
                $match_macro_name!(self, inner => inner.write_vectored(bufs))
            }

            fn write_all(&mut self, buf: &[u8]) -> ::std::io::Result<()> {
                $match_macro_name!(self, inner => inner.write_all(buf))
            }

            fn write_fmt(&mut self, fmt: ::std::fmt::Arguments<'_>) -> ::std::io::Result<()> {
                $match_macro_name!(self, inner => inner.write_fmt(fmt))
            }
        }
    }
}

/// Helper macro for [`define_whichever`], responsible for parsing and
/// delegating the trait impls.
macro_rules! impl_traits {
    ($struct_name:ident, $match_macro_name:ident, $($trait:ident),* $(,)?) => {
        $(
            $crate::utils::whichever::impl_traits!(@expand $struct_name, $match_macro_name, $trait);
        )*
    };

    (@expand $struct_name:ident, $match_macro_name:ident, Read) => {
        $crate::utils::whichever::impl_read!($struct_name, $match_macro_name);
    };

    (@expand $struct_name:ident, $match_macro_name:ident, Write) => {
        $crate::utils::whichever::impl_write!($struct_name, $match_macro_name);
    };

    (@expand $struct_name:ident, $match_macro_name:ident, $other:ident) => {
        compile_error!(concat!("Unsupported trait: ", stringify!($other)));
    };
}

/// A macro to define an enum similar to Either, but with any number of
/// variants, each containing a known (not generic) type. An invocation will
/// contain:
/// 1. A name to use for an inner macro, starting with @, typically
///    match_my_struct_name
/// 2. The enum definition, with a doc comment using attribute syntax, an
///    optional visibility specifier, and the variants
/// 3. The traits to implement, using `impl Trait {}`. Currently, we support
///    `Read` and `Write`
macro_rules! define_whichever {
    (
        @$match_macro_name:ident

        $(#[$meta:meta])*
        $vis:vis enum $struct_name:ident {
            $($variant:ident($ty:ty)),+
            $(,)?
        }
        $(impl $trait:ident for $struct_name2:ident {}),*
    ) => {
        $(#[$meta])*
        $vis enum $struct_name {
            $(
                $variant($ty),
            )+
        }

        #[macro_export]
        macro_rules! $match_macro_name {
            ($value:expr, $pattern:pat => $result:expr) => {
                match $value {
                    $(
                        $struct_name::$variant($pattern) => $result,
                    )+
                }
            };
        }

        $($crate::utils::whichever::impl_traits!{$struct_name2, $match_macro_name,  $trait})*
    };
}

pub(crate) use {define_whichever, impl_read, impl_traits, impl_write};
