# Codebook

- `valley`: Pau da Lima valley ID (1/2/4)
- `X` and `Y`: coordinates (EPSG:32724 - WGS 84 / UTM zone 24S)
- `data_type`: rat index
- `outcome`: index outcome (signs: +ve/-ve, traps: +ve/-ve, plates: number +ve)
- `offset`: number of trials (signs: 1, traps: 1 or 0.5 if closed early, plates: number remaining)
- `offset_req`: 1 if trapped closed early, 0 otherwise
- `elevation`: above sea level (m)
- `dist_trash`: 3D shortest distance to a refuse dump (m)
- `lc30_prop_veg`: Proportion of land cover which is vegetation within a 30m radius
