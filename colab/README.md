# Sharing the P-body results with the Spang lab through Google Colab

Colab runs on Google's servers and cannot reach sciCORE (accessible from the UniBas network / VPN only), so plain `rsync` to Colab is not possible. The result tables are therefore pushed from sciCORE to a **SWITCHdrive** folder with `rclone` (one-way sync, an `rsync` equivalent for cloud storage), and the Colab notebook [pbody_explorer.ipynb](pbody_explorer.ipynb) reads them from a password-protected, read-only share link of that folder.

```
sciCORE tables/…/Sanity_utils/  --rclone sync (vscode node)-->  SWITCHdrive:pbody_RNA_Spang_shared/  --WebDAV share link-->  Colab
```

## One-time setup (Zavolan lab side)

1. **App Passcode for rclone.** Sign in at <https://drive.switch.ch> with Switch edu-ID → user menu → *Settings* → *Security* → *App passcodes* → name it `scicore-rclone` → *Create new app passcode*. Note the user name and password shown (they are shown once).

2. **rclone remote** on the vscode node:

   ```bash
   module load rclone/1.68.2
   rclone config create switchdrive webdav \
       url=https://drive.switch.ch/remote.php/dav/files/<APP_PASSCODE_USER_NAME> \
       vendor=nextcloud user=<APP_PASSCODE_USER_NAME> pass=<APP_PASSCODE_PASSWORD> --obscure
   rclone mkdir switchdrive:pbody_RNA_Spang_shared
   rclone lsd switchdrive:                       # check
   ```

   The passcode is stored (obscured, not encrypted) in `~/.config/rclone/rclone.conf`; keep that file private (`chmod 600`). To revoke access, delete the passcode on drive.switch.ch.

3. **First sync**, dry run first:

   ```bash
   cd ~/Projects/pbody_RNA_Spang
   WF/scripts/sync_to_switchdrive.sh --dry-run
   WF/scripts/sync_to_switchdrive.sh
   ```

   The files that are shared are listed in [WF/scripts/sync_to_switchdrive.filter](../WF/scripts/sync_to_switchdrive.filter) (fracSanity tables, gene features, gene-set results and the tables written by `pbody_RNA_Spang.fracSanity.ipynb`; ~100 MB). `rclone sync` makes the SWITCHdrive folder identical to the selected local files, including deletions.

4. **Share link.** On drive.switch.ch, folder `pbody_RNA_Spang_shared` → *Share* → *Share link* → permissions *View only*, set a password and an expiry date. The link has the form `https://drive.switch.ch/index.php/s/<token>`. Send the link (token) and, through a separate channel, the password to the collaborators.

   Alternative: share the folder with named persons (their Switch edu-ID e-mail). Then each collaborator needs an App Passcode of their own, and the notebook's `ENDPOINTS` must use `https://drive.switch.ch/remote.php/dav/files/<their user>/pbody_RNA_Spang_shared/` with their passcode instead of the public-link endpoints.

5. **Colab notebook.** Upload [pbody_explorer.ipynb](pbody_explorer.ipynb) to Google Drive (or open it in Colab via *File → Upload notebook* and *File → Save a copy in Drive*), then *Share* it with the collaborators (*Viewer* is enough: each of them runs their own copy).

## Keeping it in sync

Run the sync after each update of the analysis: at the end of `pbody_RNA_Spang.fracSanity.ipynb`, or by hand. To sync automatically, use a cron job **on the vscode node** (compute nodes have no internet, so `scrontab` jobs would not reach SWITCHdrive; check with scicore-admin that user cron jobs are allowed there):

```bash
crontab -e
# every night at 02:30
30 2 * * * /scicore/home/zavolan/mirono0000/Projects/pbody_RNA_Spang/WF/scripts/sync_to_switchdrive.sh >> /scicore/home/zavolan/mirono0000/Projects/pbody_RNA_Spang/colab/sync.log 2>&1
```

`rclone sync --checksum` uploads only changed files, so a nightly run without changes costs a few seconds.

## Collaborator side (Spang lab)

1. Open the shared notebook in Colab and save a copy (*File → Save a copy in Drive*).
2. Add two Colab secrets (key icon in the left sidebar, switch on *Notebook access*): `PBODY_SHARE_TOKEN` (the `<token>` part of the link) and `PBODY_SHARE_PASSWORD`.
3. *Runtime → Run all*. Edit `GENES` in section 2 (gene names or systematic names) and re-run the cells below it; section 6 searches gene sets by name and lists the genes of a set.

The first run downloads the tables into the Colab session; later runs in the same session reuse them (`fetch(..., refresh=True)` forces a new download after an update).

## Notes

- The public-link WebDAV endpoint is `https://drive.switch.ch/public.php/dav/files/<token>/` (user `anonymous`, the share password) on current Nextcloud versions, `https://drive.switch.ch/public.php/webdav/` (user = token) on older ones; the notebook tries both. Not yet tested against a real SWITCHdrive share.
- Backup option if SWITCHdrive does not work out: a private Hugging Face storage bucket (`hf sync <dir> hf://buckets/<org>/<name> --delete`); data then reside in the US.
- sciCORE's own Nextcloud (io.scicore) currently refuses external access; ask scicore-admin if it becomes an option.
