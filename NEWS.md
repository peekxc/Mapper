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