#!/bin/bash

# ── Motif Scanner: Gmail SMTP setup ────────────────────────────────────────
# Run this once on each machine before using the "email results" feature.
# It writes a local .env file (not synced/copied automatically between
# machines) with the credentials app.py reads via python-dotenv.
# ──────────────────────────────────────────────────────────────────────────

cd "$(dirname "$0")"
ENV_FILE=".env"

if [ -f "$ENV_FILE" ]; then
    read -p "  .env already exists here. Overwrite? [y/N] " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo "  Cancelled."
        exit 0
    fi
fi

echo ""
echo "  Motif Scanner -- Gmail SMTP setup"
echo "  You need a Gmail App Password (not your normal Gmail password)."
echo "  Create one at: https://myaccount.google.com/apppasswords"
echo ""

read -p "  Gmail address (SMTP_USER): " smtp_user
read -s -p "  Gmail App Password (SMTP_PASS, input hidden): " smtp_pass
echo ""
read -p "  SMTP host [smtp.gmail.com]: " smtp_host
smtp_host=${smtp_host:-smtp.gmail.com}
read -p "  SMTP port [587]: " smtp_port
smtp_port=${smtp_port:-587}

cat > "$ENV_FILE" <<EOF
SMTP_HOST=$smtp_host
SMTP_PORT=$smtp_port
SMTP_USER=$smtp_user
SMTP_PASS=$smtp_pass
EOF

chmod 600 "$ENV_FILE"

echo ""
echo "  Saved to $(pwd)/$ENV_FILE (permissions locked to your user only)."
echo "  Restart the app (./launch.sh) for it to take effect."
echo ""
