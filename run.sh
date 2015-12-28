#!/bin/bash

if [ "$DEV" != "true" ]
then
  gunicorn -b "0.0.0.0:$PORT" app:app --reload
else
  python app.py
fi
