#!/bin/bash

if [ "$DEV" != "true" ]
then
  gunicorn -b "0.0.0.0:$PORT" app:app --workers=8 --threads=8 --timeout=600 --worker-class=eventlet --reload
else
  python app.py
fi
