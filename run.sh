#!/bin/bash

if [ "$DEV" != "true" ]
then
  sudo gunicorn -b "0.0.0.0:$PORT" app:app --workers=3 --timeout=600 --worker-class=eventlet --reload
else
  python app.py
fi
