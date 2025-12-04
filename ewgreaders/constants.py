Qroot = "Q:/Messdaten/Aphys_Hypothesis_data"

# metadata
md_paths = {
    "mooring": [
        "{Qroot}/{lake}/{year}/Mooring/{date}/{location}_md.json"
    ],
    "ctd": [
        "{Qroot}/{lake}/{year}/CTD/{date}/L0/{fname}.meta"
    ]
}

L0_paths = {
    "adcp": [
        "{Qroot}/{lake}/{year}/Mooring/{date}/{location}/L0/{fname}.000"
    ],
    "thermistor": [
        "{Qroot}/{lake}/{year}/Mooring/{date}/{location}/L0/{fname}.rsk"
    ],
    "o2": [
        "{Qroot}/{lake}/{year}/Mooring/{date}/{location}/L0/{fname}.rsk",
        "{Qroot}/{lake}/{year}/Mooring/{date}/{location}/L0/7450-{serial_id}/Cat.txt"           # .TXT?
    ],
    "ctd": [
        "{Qroot}/{lake}/{year}/CTD/{date}/L0/{fname}.TOB",
        "{Qroot}/{lake}/{year}/CTD/{date}/L0/{fname}.cnv",
        "{Qroot}/{lake}/{year}/CTD/{date}/L0/{fname}.rsk",
        "{Qroot}/{lake}/{year}/CTD/{date}/L0/{fname}.nc"               # split RBR profiles
    ],
    "vmp": [
        "{Qroot}/{lake}/{year}/Microstructure/{date}/L0/{fname}.p"
    ]
}

L1_paths = {
    "ctd": [
        "{Qroot}/{lake}/{year}/CTD/{date}/L1/{fname}.nc"
    ]
}