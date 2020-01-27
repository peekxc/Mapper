# Mapper 0.9.3
- Removed active X field in favor of 'data' accessor
- Added support for precomputed 'dist' objects
- Refactored covers to 'construct'; moved 'construct_cover' as higher-order method for MapperRef
- Improved BallCover to match more intuitive notion of a ball cover
- Overhauled nerve calculation to use a combination tree to speed up intersection calculations 
- Changed measure field to use proxy registry entries and by extension to accept custom distance functions 
- Added a few more unit tests

# Mapper 0.9.2
- Added use_data to MapperRef class for uniformity with other methods
- Added first tutorial vignette on shape comparison
- Fixed bug related to user other covers 

# Mapper 0.9.1
- Moved filter aspect from cover to MapperRef class, added 9 popular filter functions and a 'use_filter' API
- Can now construct arbitrary dimensional simplicial complexes
- Removed grapher visualization method, added pixiplex as a replacement
- Distinguished computation of the 'pullback' from the 'nerve' towards a cleaner modularity of components

# Mapper 0.9.0
- Added many API-breaking changes, including: 
- (1) Changed field name 'ls\_vertex\_map' to pullback to match theory more correctly
- (2) Changed name of 'rectangular'-type covers to 'interval'
- (3) Added ability to construct individual subsets of the cover via 'pullback' method
- Ported simplex tree implementation to standalone package, which is now linked externally  

# Mapper 0.8.0
- Initial release w/ preliminary code and documentation. 