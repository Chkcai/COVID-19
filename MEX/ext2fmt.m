function fmt=ext2fmt(ext)
switch ext
    case '.img'
        fmt='HFA';
    case '.tif'
        fmt='GTiff';
end
end