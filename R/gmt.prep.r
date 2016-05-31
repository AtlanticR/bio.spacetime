  # general setup of GMT output config

  gmt.prep = function(fnt="18p", psres=300) {
    cmd( "gmtset FONT_ANNOT_PRIMARY Helvetica")
      # other good: Times-Roman Palatino-Roman Helvetica-Narrow
      #             NewCenturySchlbk-Roman Courier
    cmd( "gmtset FONT_ANNOT_PRIMARY", fnt)
    cmd( "gmtset FONT_LABEL ", fnt)
#    cmd( "gmtset DOTS_PR_INCH", psres)  # obsolete in gmt5
    cmd( "gmtset PS_MEDIA letter")
    cmd( "gmtset PS_PAGE_ORIENTATION portrait")
    cmd( "gmtset PS_CHAR_ENCODING = Standard+")
#    cmd( "gmtset PSIMAGE_FORMAT bin") # binaries are smaller  obsolete in GMT5
    cmd( "gmtset PROJ_ELLIPSOID WGS-84")
    cmd( "gmtset FORMAT_GEO_MAP dddF") # add degree notation in range -180 to +180,
    cmd( "gmtset MAP_FRAME_TYPE fancy") # the + adds rounded corners
    cmd( "gmtset MAP_FRAME_WIDTH 1p") #map frame
    cmd( "gmtset MAP_FRAME_PEN 1p")
    return (NULL)
  }


