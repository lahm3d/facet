import rasterio
import rasterio.mask
import fiona
import numpy as np
import pandas as pd

def delineate(
    hand, sub_watersheds_poly, reach_id, flood_extent_layer, flood_height_thresholds, min_da, max_da, logger
):
    """
    Delineate a FIM from the HAND grid using depth at each polygon (eg, catchment)

    Args:
        hand:
        sub_watersheds_poly:
        reach_id:
        flood_extent_layer:
        flood_height_thresholds:
        logger:

    Returns:
    """
    # Open the HAND layer:
    with rasterio.open(hand) as ds_hand:

        out_meta = ds_hand.meta.copy()
        arr_fim = np.empty(
            [out_meta["height"], out_meta["width"]], dtype=out_meta["dtype"]
        )
        arr_fim[:, :] = out_meta["nodata"]

        lst_h = []
        lst_linkno = []
        lst_prov = []
        lst_da = []

        # Open the catchment polygon layer:
        with fiona.open(sub_watersheds_poly, "r") as sheds:
            for shed in sheds:
                # Get the linkno:
                linkno = shed["properties"]["LINKNO"]
                # Get the Province:
                prov = shed["properties"]["PROVINCE"]
                # Get the Drainage Area in km^2:
                da_km2 = shed["properties"]["DSContArea"] / 1000000

                lst_da.append(da_km2)

                if prov == "COASTAL PLAIN" and da_km2 >= min_da and da_km2 <= max_da:
                    h = 1.65
                elif prov == "PIEDMONT" and da_km2 >= min_da and da_km2 <= max_da:
                    h = (np.log10(da_km2) * 0.471 + 0.523) ** 2
                elif prov == "VALLEY AND RIDGE" and da_km2 >= min_da and da_km2 <= max_da:
                    h = (np.log10(da_km2) * 0.471 + 0.375) ** 2
                elif prov == "APPALACHIAN PLATEAUS" and da_km2 >= min_da and da_km2 <= max_da:
                    h = (np.log10(da_km2) * 0.471 + 0.041) ** 2
                elif prov == "BLUE RIDGE" and da_km2 >= min_da and da_km2 <= max_da:
                    # place holder hand method for blue ridge
                    # h = (np.log10(da_km2) * 0.471 + 0.041) ** 2
                    h = 1.56
                else:
                    lst_h.append(-9999)
                    lst_linkno.append(linkno)
                    lst_prov.append(prov)
                    continue  # skip this catchment

                lst_h.append(h)
                lst_linkno.append(linkno)
                lst_prov.append(prov)

                try:
                    # Mask the bankpts file for each feature:
                    w, out_transform = rasterio.mask.mask(
                        ds_hand, [shed["geometry"]], crop=True
                    )
                    w[(w > h)] = out_meta["nodata"]  # Assign NoData to everywhere else

                    # Now write out the FIM for this shed:
                    w = w[0]
                    shape = np.shape(w)

                    # window bounds in x-y space (west, south, east, north)
                    bounds = rasterio.transform.array_bounds(
                        shape[0], shape[1], out_transform
                    )

                    col_min, row_min = ~ds_hand.transform * (bounds[0], bounds[3])

                    row_min = int(row_min)
                    col_min = int(col_min)
                    row_max = int(row_min + shape[0])
                    col_max = int(col_min + shape[1])

                    arr_w = np.empty(
                        [row_max - row_min, col_max - col_min], dtype=out_meta["dtype"]
                    )
                    arr_w[:, :] = arr_fim[row_min:row_max, col_min:col_max]
                    #
                    inds_lt = np.where(arr_fim[row_min:row_max, col_min:col_max] < w)
                    arr_w[inds_lt] = w[inds_lt]

                    # assign the FIM window for this catchment to the total array
                    arr_fim[row_min:row_max, col_min:col_max] = arr_w
                except:
                    logger.info(
                        f"""
                    WARNING: Problem masking HAND grid using catchment Linkno: {linkno}
                    """
                    )

    # Write out the final FIM grid
    out_meta.update(compress="lzw")
    with rasterio.open(
        flood_extent_layer, "w", tiled=True, blockxsize=512, blockysize=512, **out_meta
    ) as dest:
        dest.write(arr_fim, indexes=1)

    # Write HAND heights to csv
    df_h = pd.DataFrame({reach_id: lst_linkno, "prov": lst_prov, "h": lst_h, "da": lst_da})
    df_h.to_csv(flood_height_thresholds)

    return